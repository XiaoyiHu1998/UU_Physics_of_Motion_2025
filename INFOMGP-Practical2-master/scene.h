#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readMESH.h>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"
#include "constraints.h"

using namespace Eigen;
using namespace std;


void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p);
void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir);
void center(const void *_obj, ccd_vec3_t *dir);



//Impulse is defined as a pair <position, direction>
typedef std::pair<RowVector3d,RowVector3d> Impulse;


//the class the contains each individual rigid objects and their functionality
class Mesh{
public:
  MatrixXd origV;   //original vertex positions, where COM=(0.0,0.0,0.0) - never change this!
  MatrixXd currV;   //current vertex position
  MatrixXi F;   //faces of the tet mesh
  MatrixXi T;   //Tets in the tet mesh
  
  VectorXi boundTets;  //indices (from T) of just the boundary tets, for collision
  
  //position of object in space. We must always have that currV = QRot(origV, orientation)+ COM
  RowVector4d orientation; //current orientation
  RowVector3d COM;  //current center of mass
  Matrix3d invIT;  //Original *inverse* inertia tensor around the COM, defined in the rest state to the object (so to the canonical world system)
  
  VectorXd tetVolumes;    //|T|x1 tetrahedra volumes
  VectorXd invMasses;     //|T|x1 tetrahedra *inverse* masses
  
  //kinematics
  bool isFixed;  //is the object immobile
  double totalMass;  //sum(1/invMass)
  double totalVolume;
  RowVector3d comVelocity;  //the linear velocity of the center of mass
  RowVector3d angVelocity;  //the angular velocity of the object.

  //checking collision between bounding boxes, and consequently the boundary tets if succeeds.
  //you do not need to update these functions (isBoxCollide and isCollide) unless you are doing a different collision
  
  bool isBoxCollide(const Mesh& m){
    RowVector3d VMin1=currV.colwise().minCoeff();
    RowVector3d VMax1=currV.colwise().maxCoeff();
    RowVector3d VMin2=m.currV.colwise().minCoeff();
    RowVector3d VMax2=m.currV.colwise().maxCoeff();
    
    //checking all axes for non-intersection of the dimensional interval
    for (int i=0;i<3;i++)
      if ((VMax1(i)<VMin2(i))||(VMax2(i)<VMin1(i)))
        return false;
    
    return true;  //all dimensional intervals are overlapping = intersection
    
  }
  
  bool isCollide(const Mesh& m, double& depth, RowVector3d& intNormal, RowVector3d& intPosition){
    
    
    if ((isFixed && m.isFixed))  //collision does nothing
      return false;
    
    //collision between bounding boxes
    if (!isBoxCollide(m))
      return false;
    
    //otherwise, full test
    ccd_t ccd;
    CCD_INIT(&ccd);
    ccd.support1       = support; // support function for first object
    ccd.support2       = support; // support function for second object
    ccd.center1         =center;
    ccd.center2         =center;
    
    ccd.first_dir       = stub_dir;
    ccd.max_iterations = 100;     // maximal number of iterations
    
    
    void* obj1=(void*)this;
    void* obj2=(void*)&m;
    
    ccd_real_t _depth;
    ccd_vec3_t dir, pos;
    
    int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);
    
    if (nonintersect)
      return false;
    
    for (int k=0;k<3;k++){
      intNormal(k)=dir.v[k];
      intPosition(k)=pos.v[k];
    }
    
    depth =_depth;
    intPosition-=depth*intNormal/2.0;
    
    //Vector3d p1=intPosition+depth*intNormal;
    //Vector3d p2=intPosition;
    //std::cout<<"intPosition: "<<intPosition<<std::endl;
    
    //std::cout<<"depth: "<<depth<<std::endl;
    //std::cout<<"After ccdGJKIntersect"<<std::endl;
    
    //return !nonintersect;
    
    return true;
    
  }
  
  
  //return the current inverted inertia tensor around the current COM. Update it by applying the orientation
  Matrix3d getCurrInvInertiaTensor(){
      Matrix3d R = Q2RotMatrix(orientation);
      return R * invIT * R.transpose();
  }
  
  
  //Update the current position and orientation by integrating the linear and angular velocities, and update currV accordingly
  //You need to modify this according to its purpose
  void updatePosition(double timeStep){
    //just forward Euler now
    if (isFixed)
      return;  //a fixed object is immobile

    //Update the center of mass position
    COM += comVelocity * timeStep; //comVelocity is alreadyt updated for t + deltaT, so therefore this is semi-implicit euler time integration

    // Calculate the angle of rotation
    double theta = angVelocity.norm() * timeStep;

    double scalarComponent = cos(theta / 2.0);
    RowVector3d vectorComponent = angVelocity.normalized() * sin(theta / 2.0);

    Quaternion<double> orientationQuaternion = Quaternion<double>(orientation[0], orientation[1], orientation[2], orientation[3]);
    Quaternion<double> rotationQuaternion = Quaternion<double>(AngleAxisd(theta, angVelocity.normalized()));

    // Update the stored orientation
    //orientationQuaternion = rotationQuaternion * orientationQuaternion * rotationQuaternion.inverse();
    orientationQuaternion = rotationQuaternion * orientationQuaternion;
    orientation[0] = orientationQuaternion.w();
    orientation[1] = orientationQuaternion.x();
    orientation[2] = orientationQuaternion.y();
    orientation[3] = orientationQuaternion.z();
    
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
  }
  
  
  RowVector3d initStaticProperties(const double density)
  {
    tetVolumes.conservativeResize(T.rows());
    
    RowVector3d naturalCOM; naturalCOM.setZero();
    Matrix3d IT; IT.setZero();
    for (int i=0;i<T.rows();i++){
      Vector3d e01=origV.row(T(i,1))-origV.row(T(i,0));
      Vector3d e02=origV.row(T(i,2))-origV.row(T(i,0));
      Vector3d e03=origV.row(T(i,3))-origV.row(T(i,0));
      Vector3d tetCentroid=(origV.row(T(i,0))+origV.row(T(i,1))+origV.row(T(i,2))+origV.row(T(i,3)))/4.0;
      tetVolumes(i)=std::abs(e01.dot(e02.cross(e03)))/6.0;
      
      naturalCOM+=tetVolumes(i)*tetCentroid;
      
    }
    
    totalVolume=tetVolumes.sum();
    totalMass=density*totalVolume;
    naturalCOM.array()/=totalVolume;
    
    //computing inertia tensor
    for (int i=0;i<T.rows();i++){
      RowVector4d xvec; xvec<<origV(T(i,0),0)-naturalCOM(0),origV(T(i,1),0)-naturalCOM(0),origV(T(i,2),0)-naturalCOM(0),origV(T(i,3),0)-naturalCOM(0);
      RowVector4d yvec; yvec<<origV(T(i,0),1)-naturalCOM(1),origV(T(i,1),1)-naturalCOM(1),origV(T(i,2),1)-naturalCOM(1),origV(T(i,3),1)-naturalCOM(1);
      RowVector4d zvec; zvec<<origV(T(i,0),2)-naturalCOM(2),origV(T(i,1),2)-naturalCOM(2),origV(T(i,2),2)-naturalCOM(2),origV(T(i,3),2)-naturalCOM(2);
      
      double I00, I11, I22, I12, I21, I01, I10, I02, I20;
      Matrix4d sumMat=Matrix4d::Constant(1.0)+Matrix4d::Identity();
      I00 = density*6*tetVolumes(i)*(yvec*sumMat*yvec.transpose()+zvec*sumMat*zvec.transpose()).sum()/120.0;
      I11 = density*6*tetVolumes(i)*(xvec*sumMat*xvec.transpose()+zvec*sumMat*zvec.transpose()).sum()/120.0;
      I22 = density*6*tetVolumes(i)*(xvec*sumMat*xvec.transpose()+yvec*sumMat*yvec.transpose()).sum()/120.0;
      I12 = I21 = -density*6*tetVolumes(i)*(yvec*sumMat*zvec.transpose()).sum()/120.0;
      I10 = I01 = -density*6*tetVolumes(i)*(xvec*sumMat*zvec.transpose()).sum()/120.0;
      I20 = I02 = -density*6*tetVolumes(i)*(xvec*sumMat*yvec.transpose()).sum()/120.0;
      
      Matrix3d currIT; currIT<<I00, I01, I02,
      I10, I11, I12,
      I20, I21, I22;
      
      IT+=currIT;
      
    }
    invIT=IT.inverse();
    if (isFixed)
      invIT.setZero();  //infinite resistance to rotation
  
    return naturalCOM;
    
  }
  
  
  //Updating the linear and angular velocities of the object
  //You need to modify this to integrate from acceleration in the field (basically gravity)
  void updateVelocity(double timeStep){
    
    if (isFixed)
      return;

    //integrating external forces (only gravity)
    Vector3d gravity; gravity << 0, -9.8, 0.0;
    comVelocity += gravity * timeStep;
  }
  
  
  //the full integration for the time step (velocity + position)
  //You need to modify this if you are changing the integration
  void integrate(double timeStep){
    updateVelocity(timeStep);
    updatePosition(timeStep);
  }
  
  
  Mesh(const MatrixXd& _V, const MatrixXi& _F, const MatrixXi& _T, const double density, const bool _isFixed, const RowVector3d& _COM, const RowVector4d& _orientation){
    origV=_V;
    F=_F;
    T=_T;
    isFixed=_isFixed;
    COM=_COM;
    orientation=_orientation;
    comVelocity.setZero();
    angVelocity.setZero();
    
    RowVector3d naturalCOM;  //by the geometry of the object
    
    //initializes the original geometric properties (COM + IT) of the object
    naturalCOM = initStaticProperties(density);
    
    origV.rowwise()-=naturalCOM;  //removing the natural COM of the OFF file (natural COM is never used again)
    
    currV.resize(origV.rows(), origV.cols());
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
    
    
    VectorXi boundVMask(origV.rows());
    boundVMask.setZero();
    for (int i=0;i<F.rows();i++)
      for (int j=0;j<3;j++)
        boundVMask(F(i,j))=1;
    
    //cout<<"boundVMask.sum(): "<<boundVMask.sum()<<endl;
    
    vector<int> boundTList;
    for (int i=0;i<T.rows();i++){
      int incidence=0;
      for (int j=0;j<4;j++)
        incidence+=boundVMask(T(i,j));
      if (incidence>2)
        boundTList.push_back(i);
    }
    
    boundTets.resize(boundTList.size());
    for (int i=0;i<boundTets.size();i++)
      boundTets(i)=boundTList[i];
  }
  
  ~Mesh(){}
};

//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
  double currTime;
  int numFullV, numFullT;
  std::vector<Mesh> meshes;
  
  //Practical 2
  vector<Constraint> constraints;   //The (user) constraints of the scene
  
  //adding an objects. You do not need to update this generally
  void addMesh(const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const double density, const bool isFixed, const RowVector3d& COM, const RowVector4d& orientation){
    
    Mesh m(V,F, T, density, isFixed, COM, orientation);
    meshes.push_back(m);
  }
  
  
  
  /*********************************************************************
   This function handles collision constraints between objects m1 and m2 when found
   Input: meshes m1, m2
   depth: the depth of penetration
   contactNormal: the normal of the conact measured m1->m2
   penPosition: a point on m2 such that if m2 <= m2 + depth*contactNormal, then penPosition+depth*contactNormal is the common contact point
   CRCoeff: the coefficient of restitution
   
   You should create a "Constraint" class, and use its resolveVelocityConstraint() and resolvePositionConstraint() *alone* to resolve the constraint.
   You are not allowed to use practical 1 collision handling
   *********************************************************************/
  void handleCollision(Mesh& m1, Mesh& m2,const double& depth, const RowVector3d& contactNormal,const RowVector3d& penPosition, const double CRCoeff, const double tolerance){
    double invMass1 = (m1.isFixed ? 0.0 : 1.0 / m1.totalMass);  //fixed meshes have infinite mass
    double invMass2 = (m2.isFixed ? 0.0 : 1.0 / m2.totalMass);

    RowVector3d contactPosition = penPosition + depth * contactNormal;

    // Positions
    Constraint collisionConstraint = Constraint(ConstraintType::COLLISION, ConstraintEqualityType::INEQUALITY, 0, m1.currV.rows(), 1, m2.currV.rows(), invMass1, invMass2, contactNormal.normalized(), -depth, CRCoeff);
    MatrixXd comPositions = Eigen::MatrixXd::Zero(2, 3);
    comPositions.row(0) << m1.COM[0], m1.COM[1], m1.COM[2];
    comPositions.row(1) << m2.COM[0], m2.COM[1], m2.COM[2];

    MatrixXd currConstPositions = MatrixXd::Zero(2, 3);
    currConstPositions.block(0, 0, 1, 3) = contactPosition;
    currConstPositions.block(1, 0, 1, 3) = penPosition;

    MatrixXd correctedPositions = Eigen::MatrixXd::Zero(2, 3);
    bool positionsAlreadyCorrect = collisionConstraint.resolvePositionConstraint(comPositions, currConstPositions, correctedPositions, tolerance);
    
    //if (!positionsAlreadyCorrect)
    {
        // Apply position resolving matrices
        m1.COM[0] = correctedPositions.row(0)[0];
        m1.COM[1] = correctedPositions.row(0)[1];
        m1.COM[2] = correctedPositions.row(0)[2];
    
        m2.COM[0] = correctedPositions.row(1)[0];
        m2.COM[1] = correctedPositions.row(1)[1];
        m2.COM[2] = correctedPositions.row(1)[2];
    }

    // Velocities
    MatrixXd currVertexPositions = MatrixXd::Zero(2, 3);
    currVertexPositions.block(0, 0, 1, 3) = contactPosition;
    currVertexPositions.block(1, 0, 1, 3) = penPosition;

    MatrixXd currComVelocities = MatrixXd::Zero(2, m1.comVelocity.cols());
    currComVelocities.row(0) = m1.comVelocity;
    currComVelocities.row(1) = m2.comVelocity;

    MatrixXd currAngularVelocities = MatrixXd::Zero(2, m1.angVelocity.cols());
    currAngularVelocities.row(0) = m1.angVelocity;
    currAngularVelocities.row(1) = m2.angVelocity;

    MatrixXd correctedCOMVelocities = MatrixXd::Zero(2, 3);
    MatrixXd correctedAngularVelocities = MatrixXd::Zero(2, 3);
    bool velocitiesAlreadyCorrect = collisionConstraint.resolveVelocityConstraint(correctedPositions, currVertexPositions, currComVelocities, 
                                                currAngularVelocities, m1.getCurrInvInertiaTensor(), m2.getCurrInvInertiaTensor(),
                                                correctedCOMVelocities, correctedAngularVelocities, tolerance);

    // Update velocity and position with corrected values
    //if (!velocitiesAlreadyCorrect)
    {
        m1.comVelocity = correctedCOMVelocities.row(0);
        m2.comVelocity = correctedCOMVelocities.row(1);
        m1.angVelocity = correctedAngularVelocities.row(0);
        m2.angVelocity = correctedAngularVelocities.row(1);
    }
  }
  
  /*********************************************************************
   This function handles a single time step by:
   1. Integrating velocities, positions, and orientations by the timeStep
   2. (Practical 2) Detecting collisions and encoding them as constraints
   3. (Practical 2) Iteratively resolved positional and velocity constraints
   
   You do not need to update this function in Practical 2
   *********************************************************************/
  void updateScene(const double timeStep, const double CRCoeff, const double tolerance, const int maxIterations){
    
    //integrating velocity, position and orientation from forces and previous states
    for (int i=0;i<meshes.size();i++)
      meshes[i].integrate(timeStep);
    

    //detecting and handling collisions when found
    //This is done exhaustively: checking every two objects in the scene.
    double depth;
    RowVector3d contactNormal, penPosition;
    for (int i=0;i<meshes.size();i++)
      for (int j=i+1;j<meshes.size();j++)
        if (meshes[i].isCollide(meshes[j],depth, contactNormal, penPosition))
          handleCollision(meshes[i], meshes[j], depth, contactNormal, penPosition, CRCoeff, tolerance);
    
    
    //Resolving user constraints iteratively until either:
    //1. Positions or velocities are valid up to tolerance (a full streak of validity in the iteration)
    //2. maxIterations has run out
    
    
    //Resolving velocity
    int currIteration=0;
    int zeroStreak=0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
    int currConstIndex=0;
    while ((zeroStreak<constraints.size())&&(currIteration*constraints.size()<maxIterations)){
      
      Constraint currConstraint=constraints[currConstIndex];
      
      RowVector3d origConstPos1=meshes[currConstraint.m2].origV.row(currConstraint.v1);
      RowVector3d origConstPos2=meshes[currConstraint.m2].origV.row(currConstraint.v2);
      
      RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
      RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;
      //cout<<"(currConstPos1-currConstPos2).norm(): "<<(currConstPos1-currConstPos2).norm()<<endl;
      //cout<<"(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm(): "<<(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm()<<endl;
      MatrixXd currCOMPositions(2,3); currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
      MatrixXd currConstPositions(2,3); currConstPositions<<currConstPos1, currConstPos2;
      MatrixXd currCOMVelocities(2,3); currCOMVelocities<<meshes[currConstraint.m1].comVelocity, meshes[currConstraint.m2].comVelocity;
      MatrixXd currAngVelocities(2,3); currAngVelocities<<meshes[currConstraint.m1].angVelocity, meshes[currConstraint.m2].angVelocity;
      
      Matrix3d invInertiaTensor1=meshes[currConstraint.m1].getCurrInvInertiaTensor();
      Matrix3d invInertiaTensor2=meshes[currConstraint.m2].getCurrInvInertiaTensor();
      MatrixXd correctedCOMVelocities, correctedAngVelocities, correctedCOMPositions;
      
      bool velocityWasValid=currConstraint.resolveVelocityConstraint(currCOMPositions, currConstPositions, currCOMVelocities, currAngVelocities, invInertiaTensor1, invInertiaTensor2, correctedCOMVelocities,correctedAngVelocities, tolerance);
      
      if (velocityWasValid){
        zeroStreak++;
      }else{
        //only update the COM and angular velocity, don't both updating all currV because it might change again during this loop!
        zeroStreak=0;
        meshes[currConstraint.m1].comVelocity =correctedCOMVelocities.row(0);
        meshes[currConstraint.m2].comVelocity =correctedCOMVelocities.row(1);
        
        meshes[currConstraint.m1].angVelocity =correctedAngVelocities.row(0);
        meshes[currConstraint.m2].angVelocity =correctedAngVelocities.row(1);
        
      }
      
      currIteration++;
      currConstIndex=(currConstIndex+1)%(constraints.size());
    }
    
    if (currIteration*constraints.size()>=maxIterations)
      cout<<"Velocity Constraint resolution reached maxIterations without resolving!"<<endl;
    
    
    //Resolving position
    currIteration=0;
    zeroStreak=0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
    currConstIndex=0;
    while ((zeroStreak<constraints.size())&&(currIteration*constraints.size()<maxIterations)){
      
      Constraint currConstraint=constraints[currConstIndex];
      
      RowVector3d origConstPos1=meshes[currConstraint.m1].origV.row(currConstraint.v1);
      RowVector3d origConstPos2=meshes[currConstraint.m2].origV.row(currConstraint.v2);
      
      RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
      RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;

      MatrixXd currCOMPositions(2,3); currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
      MatrixXd currConstPositions(2,3); currConstPositions<<currConstPos1, currConstPos2;
    
      MatrixXd correctedCOMPositions;
    
      bool positionWasValid=currConstraint.resolvePositionConstraint(currCOMPositions, currConstPositions,correctedCOMPositions, tolerance);
      
      if (positionWasValid){
        zeroStreak++;
      }else{
        //only update the COM and angular velocity, don't both updating all currV because it might change again during this loop!
        zeroStreak=0;

        meshes[currConstraint.m1].COM =correctedCOMPositions.row(0);
        meshes[currConstraint.m2].COM =correctedCOMPositions.row(1);
        
      }
      
      currIteration++;
      currConstIndex=(currConstIndex+1)%(constraints.size());
    }
    
    if (currIteration*constraints.size()>=maxIterations)
      cout<<"Position Constraint resolution reached maxIterations without resolving!"<<endl;
    
    
    //updating currV according to corrected COM
    for (int i=0;i<meshes.size();i++)
      for (int j=0;j<meshes[i].currV.rows();j++)
        meshes[i].currV.row(j)<<QRot(meshes[i].origV.row(j), meshes[i].orientation)+meshes[i].COM;
    
    currTime+=timeStep;
  }
  
  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool loadScene(const std::string dataFolder, const std::string sceneFileName, const std::string constraintFileName){
    
    ifstream sceneFileHandle, constraintFileHandle;
    sceneFileHandle.open(dataFolder+std::string("/")+sceneFileName);
    if (!sceneFileHandle.is_open())
      return false;
    int numofObjects;
    
    currTime=0;
    sceneFileHandle>>numofObjects;
    for (int i=0;i<numofObjects;i++){
      MatrixXi objT, objF;
      MatrixXd objV;
      std::string MESHFileName;
      bool isFixed;
      double youngModulus, poissonRatio, density;
      RowVector3d userCOM;
      RowVector4d userOrientation;
      sceneFileHandle>>MESHFileName>>density>>youngModulus>>poissonRatio>>isFixed>>userCOM(0)>>userCOM(1)>>userCOM(2)>>userOrientation(0)>>userOrientation(1)>>userOrientation(2)>>userOrientation(3);
      userOrientation.normalize();
      igl::readMESH(dataFolder+std::string("/")+MESHFileName,objV,objT, objF);
      
      //fixing weird orientation problem
      MatrixXi tempF(objF.rows(),3);
      tempF<<objF.col(2), objF.col(1), objF.col(0);
      objF=tempF;
      
      addMesh(objV,objF, objT,density, isFixed, userCOM, userOrientation);
      cout << "COM: " << userCOM <<endl;
      cout << "orientation: " << userOrientation <<endl;
    }
    
    //Practical 2 change
    //reading intra-mesh attachment constraints
    int numofConstraints;
    constraintFileHandle.open(dataFolder+std::string("/")+constraintFileName);
    if (!constraintFileHandle.is_open())
      return false;
    constraintFileHandle>>numofConstraints;
    for (int i=0;i<numofConstraints;i++){
      int attachM1, attachM2, attachV1, attachV2;
      constraintFileHandle>>attachM1>>attachV1>>attachM2>>attachV2;
      
      double initDist=(meshes[attachM1].currV.row(attachV1)-meshes[attachM2].currV.row(attachV2)).norm();
      //cout<<"initDist: "<<initDist<<endl;
      double invMass1 = (meshes[attachM1].isFixed ? 0.0 : 1.0/meshes[attachM1].totalMass);  //fixed meshes have infinite mass
      double invMass2 = (meshes[attachM2].isFixed ? 0.0 : 1.0/meshes[attachM2].totalMass);
      constraints.push_back(Constraint(DISTANCE, EQUALITY,attachM1, attachV1, attachM2, attachV2, invMass1,invMass2,RowVector3d::Zero(), initDist, 0.0));
      
    }
    
    return true;
  }
  
  
  Scene(){}
  ~Scene(){}
};


/*****************************Auxiliary functions for collision detection. Do not need updating********************************/

/** Support function for libccd*/
void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p)
{
  // assume that obj_t is user-defined structure that holds info about
  // object (in this case box: x, y, z, pos, quat - dimensions of box,
  // position and rotation)
  //std::cout<<"calling support"<<std::endl;
  Mesh *obj = (Mesh *)_obj;
  RowVector3d p;
  RowVector3d d;
  for (int i=0;i<3;i++)
    d(i)=_d->v[i]; //p(i)=_p->v[i];
  
  d.normalize();
  //std::cout<<"d: "<<d<<std::endl;
  
  int maxVertex=-1;
  int maxDotProd=-32767.0;
  for (int i=0;i<obj->currV.rows();i++){
    double currDotProd=d.dot(obj->currV.row(i)-obj->COM);
    if (maxDotProd < currDotProd){
      maxDotProd=currDotProd;
      //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
      maxVertex=i;
    }
    
  }
  //std::cout<<"maxVertex: "<<maxVertex<<std::endl;
  
  for (int i=0;i<3;i++)
    _p->v[i]=obj->currV(maxVertex,i);
  
  //std::cout<<"end support"<<std::endl;
}

void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir)
{
  dir->v[0]=1.0;
  dir->v[1]=0.0;
  dir->v[2]=0.0;
}

void center(const void *_obj,ccd_vec3_t *center)
{
  Mesh *obj = (Mesh *)_obj;
  for (int i=0;i<3;i++)
    center->v[i]=obj->COM(i);
}








#endif
