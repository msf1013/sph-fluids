#include <igl/opengl/glfw/Viewer.h>
#include <thread>
#include "PhysicsHook.h"
#include "FluidsHook.h"
#include <igl/unproject.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

static PhysicsHook *hook = NULL;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
        hook->run();
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    hook->reset();
}

bool drawCallback(igl::opengl::glfw::Viewer &viewer)
{
    if (!hook)
        return false;

    hook->render(viewer);
    return false;
}

bool mouseWheel(igl::opengl::glfw::Viewer& viewer, float delta)
{
    return true;
}

bool mouseClickedCallback(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    if (button != 2)
        return false;

    Eigen::Vector3f pos(viewer.down_mouse_x, viewer.core.viewport[3] - viewer.down_mouse_y, 1);
    Eigen::Matrix4f model = viewer.core.view;
    Eigen::Vector3f unproj = igl::unproject(pos, model, viewer.core.proj, viewer.core.viewport);
    Eigen::Vector4f eye = viewer.core.view.inverse() * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
    Eigen::Vector3d dir;
    for (int i = 0; i < 3; i++)
        dir[i] = unproj[i] - eye[i];
    dir.normalize();
    return hook->mouseClicked(viewer, dir);
}

bool mouseReleasedCallback(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    if (button != 2)
        return false;

    Eigen::Vector3f pos(viewer.current_mouse_x, viewer.core.viewport[3] - viewer.current_mouse_y, 1);
    Eigen::Matrix4f model = viewer.core.view;
    Eigen::Vector3f unproj = igl::unproject(pos, model, viewer.core.proj, viewer.core.viewport);
    Eigen::Vector4f eye = viewer.core.view.inverse() * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
    Eigen::Vector3d dir;
    for (int i = 0; i < 3; i++)
        dir[i] = unproj[i] - eye[i];
    dir.normalize();
    return hook->mouseReleased(viewer, dir);
}

bool keyCallback(igl::opengl::glfw::Viewer &viewer, unsigned int key, int modifiers)
{
    if (key == ' ')
    {
        toggleSimulation();
        return true;
    }
    Eigen::Vector4f look4 = viewer.core.view.inverse() * Eigen::Vector4f(0, 0, 1.0, 0.0);
    Eigen::Vector4f left4 = viewer.core.view.inverse() * Eigen::Vector4f(1.0, 0.0, 0.0, 0.0);
    Eigen::Vector3f look(look4[0], look4[1], look4[2]);
    Eigen::Vector3f left(left4[0], left4[1], left4[2]);
    if (key == 'w')
    {        
        viewer.core.camera_base_translation += look;
    }
    if (key == 's')
    {
        viewer.core.camera_base_translation -= look;
    }
    if (key == 'a')
    {
        viewer.core.camera_base_translation += left;
    }
    if (key == 'd')
    {
        viewer.core.camera_base_translation -= left;
    }
    return false;
}

bool drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Control", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Run/Pause Sim", ImVec2(-1, 0)))
        {
            toggleSimulation();
        }
        if (ImGui::Button("Reset Sim", ImVec2(-1, 0)))
        {
            resetSimulation();
        }
    }
    hook->drawGUI(menu);
    return false;
}

int main(int argc, char *argv[])
{
    igl::opengl::glfw::Viewer viewer;

    hook = new FluidsHook();
    hook->reset();

    viewer.data().set_face_based(true);
    viewer.core.is_animating = true;    

    viewer.callback_key_pressed = keyCallback;
    viewer.callback_pre_draw = drawCallback;
    viewer.callback_mouse_down = mouseClickedCallback;
    viewer.callback_mouse_up = mouseReleasedCallback;
    viewer.callback_mouse_scroll = mouseWheel;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    menu.callback_draw_viewer_menu = [&]() {drawGUI(menu); };


    viewer.launch();
}
