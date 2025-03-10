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
  
  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff):
      constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), 
      m1(_m1), v1(_v1), m2(_m2), v2(_v2), 
      invMass1(_invMass1), invMass2(_invMass2), refValue(_refValue), CRCoeff(_CRCoeff)
  {
    refVector=_refVector;
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
      if (constraintEqualityType == ConstraintEqualityType::INEQUALITY && refValue > 0)
          return true;

    //RowVectorXd constGradient(12);

    // Get info per mesh out of matrices
    RowVectorXd comPositionM1 = currCOMPositions.row(0);
    RowVectorXd comPositionM2 = currCOMPositions.row(1);

    RowVectorXd contactPosition = currVertexPositions.row(0);
    //RowVectorXd contactArmM2 = currVertexPositions.row(1);

    RowVectorXd COMVelocityM1 = currCOMVelocities.row(0);
    RowVectorXd COMVelocityM2 = currCOMVelocities.row(1);

    RowVectorXd currAngularVelocityM1 = currAngularVelocities.row(0);
    RowVectorXd currAngularVelocityM2 = currAngularVelocities.row(1);

    if (constraintType == ConstraintType::COLLISION)
    {
        // Lagrange Multiplier
        RowVector3d contactArmM1 = contactPosition - comPositionM1;
        RowVector3d contactArmM2 = contactPosition - comPositionM2;

        // Jacobian
        // if buggy, fill in per element
        // If wrong, inverse refVector (contactNormal) direction from m1->m2 to m2->m1
        VectorXd jacobian = VectorXd::Zero(12);
        jacobian.segment(0, 3) = refVector;
        jacobian.segment(3, 3) = contactArmM1.cross(refVector);
        jacobian.segment(6, 3) = -refVector;
        jacobian.segment(9, 3) = -contactArmM2.cross(refVector);

        // Velocities
        VectorXd velocityVector = VectorXd::Zero(12);
        velocityVector.segment(0, 3) = COMVelocityM1;
        velocityVector.segment(3, 3) = currAngularVelocityM1;
        velocityVector.segment(6, 3) = COMVelocityM2;
        velocityVector.segment(9, 3) = currAngularVelocityM2;

        // Stop if under tolerance
        double distance = fabs(jacobian.dot(velocityVector));
        if (distance <= tolerance)
        {
            correctedCOMVelocities = currCOMVelocities;
            correctedAngularVelocities = currAngularVelocities;
            return true;
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
    }
    else
    {
        // Handle distance constraint velocities
        return false;
    }


    
    /**************
     TODO: write velocity correction procedure:
     1. If the velocity Constraint is satisfied up to tolerate ("abs(Jv)<=tolerance"), set corrected values to original ones and return true
     
     2. Otherwise, correct linear and angular velocities as learnt in class.
     
     Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
     
     Own notes:
     express FC as FC = J * lambda
     ***************/
    
     //Stub code: remove upon implementation
     correctedCOMVelocities=currCOMVelocities;
     correctedAngularVelocities=currAngularVelocities;
     return true;
     //end of stub code
  }
  
  //projects the position unto the constraint
  //returns true if constraint was already valid with "currPositions"
  bool resolvePositionConstraint( const MatrixXd& currCOMPositions,
                                  const MatrixXd& currConstPositions,
                                  MatrixXd& correctedCOMPositions,
                                  double tolerance ) 
  {
      // Constraint is already valid.
      if (constraintEqualityType == ConstraintEqualityType::INEQUALITY && refValue > 0 || fabs(refValue) <= tolerance) {
          correctedCOMPositions = currCOMPositions;
          return true;  
      }

      RowVectorXd comPositionM1 = currCOMPositions.row(0);
      RowVectorXd comPositionM2 = currCOMPositions.row(1);

      RowVectorXd correctedPositionM1 = RowVectorXd::Zero(3);
      RowVectorXd correctedPositionM2 = RowVectorXd::Zero(3);

      double massM1 = 1.0 / invMass1;
      double massM2 = 1.0 / invMass2;

      MatrixXd invMassMatrix = MatrixXd::Zero(6, 6);
      invMassMatrix.block(0, 0, 3, 3) = invMass1 * MatrixXd::Identity(3, 3);
      invMassMatrix.block(3, 3, 3, 3) = invMass2 * MatrixXd::Identity(3, 3);

      //// See Topic 6, slides 8 & 9.
      //VectorXd jacobian = RowVectorXd::Zero(6);
      //jacobian.segment(0, 3) = (comPositionM1 - comPositionM2) / (comPositionM1 - comPositionM2).norm();
      //jacobian.segment(3, 3) = -(comPositionM1 - comPositionM2) / (comPositionM1 - comPositionM2).norm();
      //RowVectorXd jacobianTransposed = jacobian.transpose();

      if (constraintType == ConstraintType::COLLISION) // Needs to be fixed
      {
          RowVectorXd penPosition = currConstPositions.row(0);
          RowVectorXd contactPosition = currCOMPositions.row(1);

          std::cout << "extract values" << std::endl;

          RowVectorXd n = (comPositionM1 - comPositionM2).normalized();
          VectorXd jacobian = VectorXd::Zero(6);
          jacobian.segment(0, 3) =  n.transpose();
          jacobian.segment(3, 3) = -n.transpose();
          //jacobian = jacobian.transpose();
          std::cout << "jacobian constructed" << std::endl;

          VectorXd deltaPosition = VectorXd::Zero(6);
          std::cout << "old deltaPosition: " << deltaPosition << std::endl;
          std::cout << "new deltaPosition: " << (penPosition - contactPosition).transpose() << std::endl;
          deltaPosition.segment(0, 3) = (penPosition - contactPosition).transpose();
          deltaPosition.segment(3, 3) = (penPosition - contactPosition).transpose();
          std::cout << "deltaPosition constructed" << std::endl;

          double lambda = -1.0 * jacobian.dot(deltaPosition) / (jacobian.transpose() * invMassMatrix * jacobian);

          VectorXd positionCorrection = invMassMatrix * jacobian.transpose() * lambda;

          //VectorXd positionCorrection = lambda * invMassMatrix * jacobian;
          std::cout << "positionCorrection constructed" << std::endl;

		  correctedPositionM1 = comPositionM1 + positionCorrection.segment(0, 3).transpose();
		  correctedPositionM2 = comPositionM2 + positionCorrection.segment(3, 3).transpose();
          std::cout << "corrected Positions" << std::endl;
	  }
      else
      {
          RowVectorXd n = (comPositionM1 - comPositionM2).normalized();

          double displacementM1 = -(massM2 / (massM1 + massM2));
          double displacementM2 =  (massM1 / (massM1 + massM2));

          if (invMass1 == 0 && invMass2 == 0)
          {
              displacementM1 = 0.0;
              displacementM2 = 0.0;
          }
          else if (invMass1 == 0)
          {
              displacementM1 = 0.0;
              displacementM2 = 1.0;
          }
          else if (invMass2 == 0)
          {
              displacementM1 = 1.0;
              displacementM2 = 0.0;
          }

          correctedPositionM1 = displacementM1 * refValue * n + comPositionM1;
          correctedPositionM2 = displacementM2 * refValue * n + comPositionM2;
      }

      correctedCOMPositions.block(0, 0, 1, 3) = correctedPositionM1;
      correctedCOMPositions.block(1, 0, 1, 3) = correctedPositionM2;

      /**************
       TODO: write position correction procedure:
       1. If the position Constraint is satisfied up to tolerate ("abs(C(p)<=tolerance"), set corrected values to original ones and return true

       2. Otherwise, correct COM position as learnt in class. Note that since this is a linear correction, correcting COM position == correcting all positions the same offset. the currConstPositions are used to measure the constraint, and the COM values are corrected accordingly to create the effect.

       Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
       ***************/

      return false;
  }
};



#endif /* constraints_h */
