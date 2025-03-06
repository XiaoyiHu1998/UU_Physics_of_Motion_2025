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
  bool resolveVelocityConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions, 
      const MatrixXd& currCOMVelocities, const MatrixXd& currAngularVelocities, const Matrix3d& invInertiaTensor1, 
      const Matrix3d& invInertiaTensor2, MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngularVelocities, double tolerance){
    
    RowVectorXd constGradient(12);

    // Get info per mesh out of matrices
    RowVectorXd comPositionM1 = currCOMPositions.row(0);
    RowVectorXd comPositionM2 = currCOMPositions.row(1);

    RowVectorXd currVertexPositionsM1 = currVertexPositions.row(0);
    RowVectorXd currVertexPositionsM2 = currVertexPositions.row(1);

    RowVectorXd COMVelocityM1 = currCOMVelocities.row(0);
    RowVectorXd COMVelocityM2 = currCOMVelocities.row(1);

    RowVectorXd currAngularVelocityM1 = currAngularVelocities.row(0);
    RowVectorXd currAngularVelocityM2 = currAngularVelocities.row(1);

    if (constraintType == ConstraintType::COLLISION)
    {
        std::cout << "\n\n\n__________________NEW CONSTRAINT RESOLVE VELOCITIES_______________" << std::endl;
        std::cout << "refValue " << refValue << std::endl;

        RowVectorXd penPosition = constraintType == ConstraintType::COLLISION ? currCOMPositions.row(2) : currCOMPositions.row(1);
        std::cout << "penPosition " << penPosition << std::endl;

        // Lagrange Multiplier
        RowVector3d contactPosition = penPosition - refValue * refVector; // refValue = -depth & refVector = contactNormal
        std::cout << "contactPosition " << contactPosition << std::endl;

        RowVector3d contactArmM1 = contactPosition - comPositionM1;
        RowVector3d contactArmM2 = contactPosition - comPositionM2;
        std::cout << "contactArmM1: " << contactArmM1 << std::endl;
        std::cout << "contactArmM2: " << contactArmM2 << std::endl;
        double enumerator = (COMVelocityM1 - COMVelocityM2).dot(refVector);

        RowVector3d termM1Factor1 = contactArmM1.cross(refVector) * invInertiaTensor1;
        Vector3d termM1Factor2 = contactArmM1.cross(refVector).transpose();
        RowVector3d termM2Factor1 = contactArmM2.cross(refVector) * invInertiaTensor2;
        Vector3d termM2Factor2 = contactArmM2.cross(refVector).transpose();

        // Explicit cross product and coefficient-wise product
        double termM1 = termM1Factor1 * termM1Factor2;
        double termM2 = termM2Factor1 * termM2Factor2;

        double denomerator = (invMass1 + invMass2) + termM1 + termM2;
        double lagrangeMultiplier = - enumerator / denomerator;
        std::cout << "lagrangeMultiplier " << lagrangeMultiplier << std::endl;

        // Jacobian
        // if buggy, fill in per element
        // If wrong, inverse refVector (contactNormal) direction from m1->m2 to m2->m1
        //refVector = -refVector;
        VectorXd jacobian = RowVectorXd::Zero(12);
        jacobian.segment(0, 3) = refVector;
        jacobian.segment(3, 3) = contactArmM1.cross(refVector);
        jacobian.segment(6, 3) = -refVector;
        jacobian.segment(9, 3) = -contactArmM2.cross(refVector);
        RowVectorXd jacobianTransposed = jacobian.transpose();
        std::cout << "jacobian " << jacobian << std::endl;

        // Velocities
        VectorXd velocityVector = VectorXd::Zero(12);
        velocityVector.segment(0, 3) = COMVelocityM1;
        velocityVector.segment(3, 3) = currAngularVelocityM1;
        velocityVector.segment(6, 3) = COMVelocityM2;
        velocityVector.segment(9, 3) = currAngularVelocityM2;

        // Stop if under tolerance
        double distance = abs(jacobian.dot(velocityVector));
        std::cout << distance << std::endl;
        if (distance <= tolerance)
        {
            correctedCOMVelocities = currCOMVelocities;
            correctedAngularVelocities = currAngularVelocities;
            std::cout << "Under tolerance" << std::endl;
            return true;
        }
        std::cout << "Over tolerance" << std::endl;

        // Inv Mass Matrix
        MatrixXd massMatrix = MatrixXd::Zero(12, 12);
        massMatrix.block(0, 0, 3, 3) = invMass1 * MatrixXd::Identity(3, 3);
        massMatrix.block(3, 3, 3, 3) = invInertiaTensor1;
        massMatrix.block(6, 6, 3, 3) = invMass2 * MatrixXd::Identity(3, 3);
        massMatrix.block(9, 9, 3, 3) = invInertiaTensor2;

        std::cout << "massMatrix\n " << massMatrix << "\n_______________________________________" << std::endl;
        
        RowVectorXd deltaVelocity = lagrangeMultiplier * massMatrix * jacobian;

        std::cout << "deltaVelocity " << deltaVelocity << std::endl;

        correctedCOMVelocities.row(0) = COMVelocityM1 + deltaVelocity.block(0, 0, 1, 3);
        correctedAngularVelocities.row(0) = currAngularVelocityM1 + deltaVelocity.block(0, 3, 1, 3);

        correctedCOMVelocities.row(1) = COMVelocityM2 + deltaVelocity.block(0, 6, 1, 3);
        correctedAngularVelocities.row(1) = currAngularVelocityM2 + deltaVelocity.block(0, 9, 1, 3);

        return false;
    }
    else
    {

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
  bool resolvePositionConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions, MatrixXd& correctedCOMPositions, double tolerance){
    
    MatrixXd invMassMatrix=MatrixXd::Zero(6,6);
    
    
    /**************
     TODO: write position correction procedure:
     1. If the position Constraint is satisfied up to tolerate ("abs(C(p)<=tolerance"), set corrected values to original ones and return true
     
     2. Otherwise, correct COM position as learnt in class. Note that since this is a linear correction, correcting COM position == correcting all positions the same offset. the currConstPositions are used to measure the constraint, and the COM values are corrected accordingly to create the effect.
     
     Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.

     Own notes:
     For collision handling: currCOMPositions (ComM1,ComM2) as 2x3 matrix, penposition, correctedCOMPositions, tolerance
     ***************/
    
    // No constraint violation
    if (abs(refValue) <= tolerance)
    {
        correctedCOMPositions = currCOMPositions;
        return true;
    }

    // Constraint violated
    if (constraintType == ConstraintType::COLLISION)
    {
        
    }

    return false;
  }
};



#endif /* constraints_h */
