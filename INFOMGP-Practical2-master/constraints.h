#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType{DISTANCE, COLLISION} ConstraintType;   //You can expand it for more constraints
typedef enum ConstraintEqualityType{EQUALITY, INEQUALITY} ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint{
public:
  
  int m1, m2;                     //Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
  int v1, v2;                     //Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
  double invMass1, invMass2;      //inverse masses of two bodies
  double refValue;                //Reference values to use in the constraint, when needed (like distance)
  RowVector3d refVector;          //Reference vector when needed (like vector)
  double CRCoeff;                 //extra velocity bias
  ConstraintType constraintType;  //The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
  ConstraintEqualityType constraintEqualityType;  //whether the constraint is an equality or an inequality
  
  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff) :constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), m1(_m1), v1(_v1), m2(_m2), v2(_v2), invMass1(_invMass1), invMass2(_invMass2), refValue(_refValue), CRCoeff(_CRCoeff) {
      refVector = _refVector;
  }
  
  ~Constraint(){}
  
  
  
  //computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
  //The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
  //returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction done)
  //currCOMPositions is a 2x3 matrix, where each row is per one of the sides of the constraints; the rest of the relevant variables are similar, and so should the outputs be resized.
  bool resolveVelocityConstraint(const MatrixXd& currCOMPositions, 
                                 const MatrixXd& currVertexPositions, 
                                 const MatrixXd& currCOMVelocities, 
                                 const MatrixXd& currAngularVelocities, 
                                 const Matrix3d& invInertiaTensor1, 
                                 const Matrix3d& invInertiaTensor2, 
                                 MatrixXd& correctedCOMVelocities, 
                                 MatrixXd& correctedAngularVelocities, 
                                 double tolerance)
  {
    // Set corrected value just in case
    correctedCOMVelocities = currCOMVelocities;
    correctedAngularVelocities = currAngularVelocities;

    // Get info per mesh out of matrices
    RowVectorXd comPositionM1 = currCOMPositions.row(0);
    RowVectorXd comPositionM2 = currCOMPositions.row(1);
    RowVectorXd COMVelocityM1 = currCOMVelocities.row(0);
    RowVectorXd COMVelocityM2 = currCOMVelocities.row(1);
    RowVectorXd currAngularVelocityM1 = currAngularVelocities.row(0);
    RowVectorXd currAngularVelocityM2 = currAngularVelocities.row(1);

    // Lagrange Multiplier
    RowVectorXd contactPosition = currVertexPositions.row(0);
    RowVectorXd penPosition = currVertexPositions.row(1);
    RowVector3d contactArmM1 = contactPosition - comPositionM1;
    RowVector3d contactArmM2 = contactPosition - comPositionM2;

    RowVector3d chainNormal = (penPosition - contactPosition);
    Vector3d referenceVector = constraintType == ConstraintType::COLLISION ? refVector : chainNormal;

    // Jacobian
    // Gradient either along the chain or penetration points depending on the case
    VectorXd jacobian = VectorXd::Zero(12);
    jacobian.segment(0, 3) = referenceVector;
    jacobian.segment(3, 3) = contactArmM1.cross(referenceVector);
    jacobian.segment(6, 3) = -referenceVector;
    jacobian.segment(9, 3) = -contactArmM2.cross(referenceVector);

    // Velocities
    VectorXd velocityVector = VectorXd::Zero(12);
    velocityVector.segment(0, 3) = COMVelocityM1;
    velocityVector.segment(3, 3) = currAngularVelocityM1;
    velocityVector.segment(6, 3) = COMVelocityM2;
    velocityVector.segment(9, 3) = currAngularVelocityM2;

    // Check if constraint needs to be applied
    if (constraintEqualityType == ConstraintEqualityType::EQUALITY && fabs(jacobian.dot(velocityVector)) <= tolerance)
    {
        correctedCOMVelocities = currCOMVelocities;
        correctedAngularVelocities = currAngularVelocities;
        return true;
    }
    else // Inequality case
    {
        // Does this even exist, Matthijs?
    }

    // Inv Mass Matrix
    MatrixXd invMassMatrix = MatrixXd::Zero(12, 12);
    invMassMatrix.block(0, 0, 3, 3) = invMass1 * MatrixXd::Identity(3, 3);
    invMassMatrix.block(3, 3, 3, 3) = invInertiaTensor1;
    invMassMatrix.block(6, 6, 3, 3) = invMass2 * MatrixXd::Identity(3, 3);
    invMassMatrix.block(9, 9, 3, 3) = invInertiaTensor2;

    double enumerator = -(1.0 + CRCoeff) * jacobian.dot(velocityVector);
    double denominator = jacobian.transpose() * invMassMatrix * jacobian;
    VectorXd deltaVelocity = (enumerator / denominator) * (invMassMatrix * jacobian);

    correctedCOMVelocities.row(0) = COMVelocityM1 + deltaVelocity.segment(0, 3).transpose();
    correctedAngularVelocities.row(0) = currAngularVelocityM1 + deltaVelocity.segment(3, 3).transpose();

    correctedCOMVelocities.row(1) = COMVelocityM2 + deltaVelocity.segment(6, 3).transpose();
    correctedAngularVelocities.row(1) = currAngularVelocityM2 + deltaVelocity.segment(9, 3).transpose();

    return false;
    
    /**************
     TODO: write velocity correction procedure:
     1. If the velocity Constraint is satisfied up to tolerate ("abs(Jv)<=tolerance"), set corrected values to original ones and return true
     
     2. Otherwise, correct linear and angular velocities as learnt in class.
     
     Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
     
     Own notes:
     express FC as FC = J * lambda
     ***************/
    
     //Stub code: remove upon implementation
     //correctedCOMVelocities=currCOMVelocities;
     //correctedAngularVelocities=currAngularVelocities;
     //return true;
     //end of stub code
  }
  
  //projects the position unto the constraint
  //returns true if constraint was already valid with "currPositions"
  bool resolvePositionConstraint( const MatrixXd& currCOMPositions,
                                  const MatrixXd& currConstPositions,
                                  MatrixXd& correctedCOMPositions,
                                  double tolerance ) 
  {
      // Inequality Case if constraint is already valid.
      if (constraintEqualityType == ConstraintEqualityType::INEQUALITY && refValue >= 0) { // Check if refValue = -depth affects this check
          correctedCOMPositions = currCOMPositions;
          return true;  
      }

      // Set positions
      RowVectorXd comPositionM1 = currCOMPositions.row(0);
      RowVectorXd comPositionM2 = currCOMPositions.row(1);
      RowVectorXd correctedPositionM1 = comPositionM1;
      RowVectorXd correctedPositionM2 = comPositionM2;
      correctedCOMPositions = currCOMPositions; // initialize values to prevent crashes on loading constraint scene

      // Inverse Mass Matrix
      double massM1 = 1.0 / invMass1;
      double massM2 = 1.0 / invMass2;

      MatrixXd invMassMatrix = MatrixXd::Zero(6, 6);
      invMassMatrix.block(0, 0, 3, 3) = invMass1 * MatrixXd::Identity(3, 3);
      invMassMatrix.block(3, 3, 3, 3) = invMass2 * MatrixXd::Identity(3, 3);

      if (constraintType == ConstraintType::COLLISION) // Needs to be fixed
      {
          RowVectorXd penPosition = currConstPositions.row(0);
          RowVectorXd contactPosition = currConstPositions.row(1);

          RowVectorXd n = (penPosition - contactPosition).normalized();
          VectorXd jacobian = VectorXd::Zero(6);
          jacobian.segment(0, 3) =  n.transpose();
          jacobian.segment(3, 3) = -n.transpose();

          VectorXd deltaPosition = VectorXd::Zero(6);
          deltaPosition.segment(0, 3) = (penPosition - contactPosition).transpose();
          deltaPosition.segment(3, 3) = (penPosition - contactPosition).transpose();

          double lambda = -1.0 * jacobian.dot(deltaPosition) / (jacobian.transpose() * invMassMatrix * jacobian);
          VectorXd positionCorrection = invMassMatrix * jacobian.transpose() * lambda;

		  correctedPositionM1 = comPositionM1 + positionCorrection.segment(0, 3).transpose();
		  correctedPositionM2 = comPositionM2 + positionCorrection.segment(3, 3).transpose();
	  }
      else // constraintType != ConstraintType::COLLISION
      {
          // currConstPositions should be penetration points in distance constraints
          double currDistance = (currConstPositions.row(0) - currConstPositions.row(1)).norm();
          double constraintValue = currDistance - refValue;

          if (fabs(constraintValue) <= tolerance) // correctedPositions already set to currCOMPositions before if block
              return true;

          RowVectorXd n = (currConstPositions.row(0) - currConstPositions.row(1)).normalized();

          double displacementM1 = -(massM2 / (massM1 + massM2));
          double displacementM2 =  (massM1 / (massM1 + massM2));

          correctedPositionM1 = displacementM1 * constraintValue * n + comPositionM1;
          correctedPositionM2 = displacementM2 * constraintValue * n + comPositionM2;
      }

      correctedCOMPositions.block(0, 0, 1, 3) = correctedPositionM1.segment(0, 3);
      correctedCOMPositions.block(1, 0, 1, 3) = correctedPositionM2.segment(0, 3);
      return false;

      /**************
       TODO: write position correction procedure:
       1. If the position Constraint is satisfied up to tolerate ("abs(C(p)<=tolerance"), set corrected values to original ones and return true

       2. Otherwise, correct COM position as learnt in class. Note that since this is a linear correction, correcting COM position == correcting all positions the same offset. the currConstPositions are used to measure the constraint, and the COM values are corrected accordingly to create the effect.

       Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
       ***************/
  }
};



#endif /* constraints_h */
