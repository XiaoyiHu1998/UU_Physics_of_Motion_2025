#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <iostream>
#include "scene.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F;
igl::opengl::glfw::Viewer mgpViewer;

float currTime = 0;
bool animationHack;  //fixing the weird camera bug in libigl
//initial values
float timeStep = 0.02;
float CRCoeff= 1.0;
double dragCoefficient = 2;
double kineticFrictionCoefficient = 0.5;
float dragCoefficient = 0.0;

Scene scene;

Eigen::MatrixXd platV;
Eigen::MatrixXi platF;
Eigen::MatrixXi platT;
Eigen::RowVector3d platCOM;
Eigen::RowVector4d platOrientation;

void createPlatform()
{
  double platWidth=100.0;
  platCOM<<0.0,-5.0,-0.0;
  platV.resize(9,3);
  platF.resize(12,3);
  platT.resize(12,4);
  platV<<-platWidth,0.0,-platWidth,
  -platWidth,0.0,platWidth,
  platWidth,0.0,platWidth,
  platWidth,0.0, -platWidth,
  -platWidth,-platWidth/10.0,-platWidth,
  -platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0, -platWidth,
  0.0,-platWidth/20.0, 0.0;
  platF<<0,1,2,
  2,3,0,
  6,5,4,
  4,7,6,
  1,0,5,
  0,4,5,
  2,1,6,
  1,5,6,
  3,2,7,
  2,6,7,
  0,3,4,
  3,7,4;
  
  platOrientation<<1.0,0.0,0.0,0.0;
  
  platT<<platF, VectorXi::Constant(12,8);
  
  
}

void updateMeshes(igl::opengl::glfw::Viewer &viewer)
{
  RowVector3d platColor; platColor<<0.8,0.8,0.8;
  RowVector3d meshColor; meshColor<<0.8,0.2,0.2;
  
  for (int i=0;i<scene.meshes.size();i++){
    viewer.data_list[i].clear();
    viewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
    viewer.data_list[i].set_face_based(true);
    viewer.data_list[i].set_colors(meshColor);
    viewer.data_list[i].show_lines=false;
  }
  viewer.data_list[0].show_lines=false;
  viewer.data_list[0].set_colors(platColor.replicate(scene.meshes[0].F.rows(),1));
  viewer.data_list[0].set_face_based(true);
  
  viewer.core().align_camera_center(scene.meshes[0].currV);
}



bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  if (key == ' ')
  {
    viewer.core().is_animating = !viewer.core().is_animating;
    if (viewer.core().is_animating)
      cout<<"Simulation running"<<endl;
    else
      cout<<"Simulation paused"<<endl;
    return true;
  }
  
  if (key == 'S')
  {
    if (!viewer.core().is_animating){
      scene.updateScene(timeStep, CRCoeff, dragCoefficient, kineticFrictionCoefficient);
      currTime+=timeStep;
      updateMeshes(viewer);
      std::cout <<"currTime: "<<currTime<<std::endl;
      return true;
    }
  }
  return false;
}




bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
  using namespace Eigen;
  using namespace std;
  
  if (viewer.core().is_animating){
    if (!animationHack)
      scene.updateScene(timeStep, CRCoeff, dragCoefficient, kineticFrictionCoefficient);
    else
      viewer.core().is_animating=false;
    animationHack=false;
    currTime+=timeStep;
    updateMeshes(viewer);
  }
 
  
  return false;
}

class CustomMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{
  
  virtual void draw_viewer_menu() override
  {
    // Draw parent menu
    //ImGuiMenu::draw_viewer_menu();
    
    // Add new group
    if (ImGui::CollapsingHeader("Algorithm Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
      ImGui::InputFloat("CR Coeff", &CRCoeff, 0, 0, "%.2f");
      if (ImGui::InputFloat("Time Step", &timeStep)) {
        mgpViewer.core().animation_max_fps = (((int)1.0/timeStep));
      }

      ImGui::InputFloat("dragCoefficient", &dragCoefficient);
    }

    if (ImGui::CollapsingHeader("Mesh Velocities", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // Draw Mesh velocities
        for (int i = 0; i < scene.meshes.size(); i++)
        {
            std::string label = "Mesh " + std::to_string(i);
            ImGui::Text(label.c_str());

            // Velocity Editing
            float velocityX = static_cast<float>(scene.meshes[i].comVelocity.data()[0]);
            float velocityY = static_cast<float>(scene.meshes[i].comVelocity.data()[1]);
            float velocityZ = static_cast<float>(scene.meshes[i].comVelocity.data()[2]);
            float velocityFloat[3] = {velocityX, velocityY, velocityZ};

            std::string velocityLabel = "M" + std::to_string(i) + " velocity";
            ImGui::DragFloat3(velocityLabel.c_str(), velocityFloat);
            scene.meshes[i].comVelocity.data()[0] = static_cast<double>(velocityFloat[0]);
            scene.meshes[i].comVelocity.data()[1] = static_cast<double>(velocityFloat[1]);
            scene.meshes[i].comVelocity.data()[2] = static_cast<double>(velocityFloat[2]);

            // Impulse Editing
            float ImpulsePositionX = static_cast<float>(scene.meshes[i].ExternalImpulsePosition.data()[0]);
            float ImpulsePositionY = static_cast<float>(scene.meshes[i].ExternalImpulsePosition.data()[1]);
            float ImpulsePositionZ = static_cast<float>(scene.meshes[i].ExternalImpulsePosition.data()[2]);
            float ImpulsePositionFloat[3] = { ImpulsePositionX, ImpulsePositionY, ImpulsePositionZ };

            std::string impulsePositionLabel = "M" + std::to_string(i) + " impulsePosition";
            if (ImGui::DragFloat3(impulsePositionLabel.c_str(), ImpulsePositionFloat))
            {
                scene.meshes[i].ExternalImpulsePosition  = RowVector3d(static_cast<double>(ImpulsePositionFloat[0]), static_cast<double>(ImpulsePositionFloat[1]), static_cast<double>(ImpulsePositionFloat[2]));
            }

            float ImpulseDirectionX = static_cast<float>(scene.meshes[i].ExternalImpulseDirection.data()[0]);
            float ImpulseDirectionY = static_cast<float>(scene.meshes[i].ExternalImpulseDirection.data()[1]);
            float ImpulseDirectionZ = static_cast<float>(scene.meshes[i].ExternalImpulseDirection.data()[2]);
            float ImpulseDirectionFloat[3] = { ImpulseDirectionX, ImpulseDirectionY, ImpulseDirectionZ };

            std::string impulseDirectionLabel = "M" + std::to_string(i) + " impulseDirection";
            if (ImGui::DragFloat3(impulseDirectionLabel.c_str(), ImpulseDirectionFloat))
            {
                scene.meshes[i].ExternalImpulseDirection = RowVector3d(static_cast<double>(ImpulseDirectionFloat[0]), static_cast<double>(ImpulseDirectionFloat[1]), static_cast<double>(ImpulseDirectionFloat[2]));
            }

            std::string impulseApplyLabel = "M" + std::to_string(i) + " Apply Impulse";
            if (ImGui::Button(impulseApplyLabel.c_str()))
            {
                RowVector3d impulsePosition = scene.meshes[i].ExternalImpulsePosition;
                RowVector3d impulseDirection = scene.meshes[i].ExternalImpulseDirection;

                scene.meshes[i].currImpulses.push_back(Impulse(impulsePosition, impulseDirection));
                scene.meshes[i].updateImpulseVelocities();

                scene.meshes[i].ExternalImpulsePosition.setZero();
                scene.meshes[i].ExternalImpulseDirection.setZero();
            }

            ImGui::Separator();
        }
      ImGui::InputDouble("dragCoefficient", &dragCoefficient);
      ImGui::InputDouble("kineticFrictionCoefficient", &kineticFrictionCoefficient);
    }
  }
};



int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  
  // Load scene
  if (argc<3){
    cout<<"Please provide path (argument 1 and name of scene file (argument 2)!"<<endl;
    return 0;
  }
  cout<<"scene file: "<<std::string(argv[2])<<endl;
  //create platform
  createPlatform();
  scene.addMesh(platV, platF, platT, 10000.0, true, platCOM, platOrientation);
  
  //load scene from file
  scene.loadScene(std::string(argv[1]),std::string(argv[2]));

  scene.updateScene(0.0, CRCoeff, dragCoefficient, kineticFrictionCoefficient);
  
  // Viewer Settings
  for (int i=0;i<scene.meshes.size();i++){
    if (i!=0)
      mgpViewer.append_mesh();
    //mgpViewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
  }
  //mgpViewer.core.align_camera_center(scene.meshes[0].currV);
  
  
  mgpViewer.callback_pre_draw = &pre_draw;
  mgpViewer.callback_key_down = &key_down;
  mgpViewer.core().is_animating = true;
  animationHack = true;
  mgpViewer.core().animation_max_fps = 50.;

  CustomMenu menu;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  mgpViewer.plugins.push_back(&plugin);
  plugin.widgets.push_back(&menu);
  
  cout<<"Press [space] to toggle continuous simulation" << endl;
  cout<<"Press 'S' to advance time step-by-step"<<endl;
  updateMeshes(mgpViewer);
  mgpViewer.launch();
  
}
