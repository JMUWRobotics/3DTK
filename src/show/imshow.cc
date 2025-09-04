/*
 * imshow implementation
 *
 * Copyright (C) Fabian Arzberger.
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief This is show, but without GLUI. Instead, ImGui is used.
 * @author Fabian Arzberger. Institute of CS, University of Wuerzburg, Germany.
 *
 * Original authors of show.cc and show_gl.cc:
 * @author Andreas Nuechter. Institute of CS, University of Osnabrueck, Germany.
 * @author Kai Lingemann. Institute of CS, University of Osnabrueck, Germany.
 * @author Jan Elseberg. Jacobs University bremen gGmbH, Germany.
 * @author Dorit Borrmann. Jacobs University bremen gGmbH, Germany.
 */

#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui.h"
#include "imgui_impl_glut.h"
#include "imgui_impl_opengl2.h"
#define GL_SILENCE_DEPRECATION
#include "imGuIZMOquat.h"

#include "show/callbacks_glut.h"
#include "show/program_options.h"
#include "show/show_common.h"
#include "show/show_gl.h"

#include <thread>
#ifdef __linux__
#include <linux/input.h>
#endif
#ifndef __APPLE__
#ifdef SPACEMOUSE
#include "show/show_gl.h"
#include <fcntl.h>
#include <linux/input.h>
#include <math.h>
#include <spnav.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>

#define SPNAV_DEAD_THRESH_R 50 // Adjust to your liking
#define SPNAV_DEAD_THRESH_T 30 // Adjust to your liking
#endif

// Spacemouse handler:
int spacenavHandlerIm()
{
  int fixRoation = 0;
  int fixTranslation = 0;
  float translationMultiplier = 0.05; // adjust to your liking
  float rotationMultipler = 0.005;
#ifdef SPACEMOUSE
  spnav_event sev;
  if (spnav_open() == -1) {
    fprintf(stderr, "failed to connect to the space navigator daemon\n");
    return 1;
  }
  // spnav_sensitivity(0.03); // Saturates any input on my spacemouse...
  // use multipliers above instead - Fabi
  while (true) {
    int event = spnav_poll_event(&sev);
    if (event != 0) {
      // Use for Debug:
      // printf("got motion event: t(%d, %d, %d) ",
      // sev.motion.x, sev.motion.y, sev.motion.z);
      // printf("r(%d, %d, %d)\n", sev.motion.rx,
      // sev.motion.ry, sev.motion.rz);
      if (sev.type == SPNAV_EVENT_MOTION) {

        // Set rotation deadzone
        if (fabs(sev.motion.rx) < SPNAV_DEAD_THRESH_R)
          sev.motion.rx = 0;
        if (fabs(sev.motion.ry) < SPNAV_DEAD_THRESH_R)
          sev.motion.ry = 0;
        if (fabs(sev.motion.rz) < SPNAV_DEAD_THRESH_R)
          sev.motion.rz = 0;
        if (fabs(sev.motion.x) < SPNAV_DEAD_THRESH_T)
          sev.motion.x = 0;
        if (fabs(sev.motion.y) < SPNAV_DEAD_THRESH_T)
          sev.motion.y = 0;
        if (fabs(sev.motion.z) < SPNAV_DEAD_THRESH_T)
          sev.motion.z = 0;

        // Handle the event
        if (fixRoation) {
          moveCamera(-sev.motion.x * translationMultiplier,
               -sev.motion.y * translationMultiplier,
               sev.motion.z * translationMultiplier, 0, 0, 0, true);
        } else if (fixTranslation) {
          moveCamera(0, 0, 0, -0.5 * sev.motion.rx * rotationMultipler,
               0.5 * sev.motion.ry * rotationMultipler,
               -0.5 * sev.motion.rz * rotationMultipler, true);
        } else
          moveCamera(
              -sev.motion.x * translationMultiplier,
              -sev.motion.y * translationMultiplier, sev.motion.z * translationMultiplier,
              -rotationMultipler * sev.motion.rx, rotationMultipler * sev.motion.ry,
              -rotationMultipler * sev.motion.rz, true);
      } else { /* SPNAV_EVENT_BUTTON */
        if (sev.button.press && sev.button.bnum == 0) {
          fixRoation = !fixRoation;
        }
        if (sev.button.press && sev.button.bnum == 1) {
          fixTranslation = !fixTranslation;
        }
      }
    }
  }
#else
#ifdef __linux__
  struct input_id device_info;
  int fd = open("/dev/input/by-id/usb-3Dconnexion_SpaceMouse_Compact-event-if00", O_RDONLY);
  if (fd == -1) {
    perror("open");
    fprintf(stderr, "3D Mouse not connected\n");
    return 1;
  }
  if (ioctl(fd, EVIOCGID, &device_info)) {
    perror("3D Mouse EVIOCGID ioctl failed");
    return 1;
  }
  if (device_info.vendor != 0x256f) {
    fprintf(stderr, "unexpected vendor: %d\n", device_info.vendor);
    return 1;
  }
  if (device_info.product != 0xc635) {
    fprintf(stderr, "unexpected product: %d\n", device_info.product);
    return 1;
  }
  if (device_info.bustype != BUS_USB) {
    fprintf(stderr, "unexpected bus type: %d\n", device_info.bustype);
    return 1;
  }
  int axes[6] = {0, 0, 0, 0, 0, 0};
  int buttons[2] = {0, 0};
  struct input_event ev;
  uint8_t evtype_bitmask[EV_MAX / 8 + 1];
  int ev_type;

  memset(evtype_bitmask, 0, sizeof(evtype_bitmask));
  if (ioctl(fd, EVIOCGBIT(0, sizeof(evtype_bitmask)), evtype_bitmask) < 0) {
    perror("3D Mouse EVIOCGBIT ioctl failed");
    return 1;
  }
  for (;;) {
    int n = read(fd, &ev, sizeof(struct input_event));
    if (n < sizeof(struct input_event)) {
      fprintf(stderr, "unexpected event size\n");
      return 1;
    }
    switch (ev.type) {
    case EV_KEY:
      if (ev.code != 0 and ev.code != 1) {
        fprintf(stderr, "unexpected button code: %d\n", ev.code);
      }
      buttons[ev.code] = ev.value;
      break;
    case EV_REL:
      switch (ev.code) {
      case 0:
        moveCamera(-ev.value / 4, 0, 0, 0, 0, 0);
        break;
      case 1:
        moveCamera(0, 0, -ev.value / 4, 0, 0, 0);
        break;
      case 2:
        moveCamera(0, ev.value / 4, 0, 0, 0, 0);
        break;
      case 3:
        moveCamera(0, 0, 0, -ev.value / 40, 0, 0);
        break;
      case 4:
        moveCamera(0, 0, 0, 0, 0, ev.value / 40);
        break;
      case 5:
        moveCamera(0, 0, 0, 0, -ev.value / 40, 0);
        break;
      }
    default:
      cout << "Neither key nor rel" << endl;
    }
  }
#endif
#endif
  return 0;
}
#endif

// Static vars to determine starting dimensions and position of the window
static int SCREEN_WIDTH;
static int SCREEN_HEIGHT;
static int START_WIDTH_IMGUI;
static int START_HEIGHT_IMGUI;
static int START_X_IMGUI;
static int START_Y_IMGUI;

// show needs these to be defined, altough being useless when using imgui
// instead of glui...
void updateCamControls() { return; }
void resetRotationButton() { return; }
void updateViewModeControls() { return; }
void updateControls() { return; }
void updatePointModeControls() { return; }
void DrawTypeLegend();
void setup_camera();
void setup_fog();

void mapColorToValueIm()
{
  switch (listboxColorVal) {
  case 0:
    cm->setCurrentType(PointType::USE_HEIGHT);
    break;
  case 1:
    cm->setCurrentType(PointType::USE_REFLECTANCE);
    break;
  case 2:
    cm->setCurrentType(PointType::USE_TEMPERATURE);
    break;
  case 3:
    cm->setCurrentType(PointType::USE_AMPLITUDE);
    break;
  case 4:
    cm->setCurrentType(PointType::USE_DEVIATION);
    break;
  case 5:
    cm->setCurrentType(PointType::USE_TYPE);
    break;
  case 6:
    cm->setCurrentType(PointType::USE_COLOR);
    break;
  default:
    break;
  };
}

// Forward declaration
void DrawPointsIm(GLenum, bool);

extern bool classLabels;

// Show Interrupt handlers:
static bool interrupted = false;

void interruptDrawing() { interrupted = true; }

void checkForInterrupt() { interrupted = false; }

bool isInterrupted()
{
#ifndef __APPLE__
  glutMainLoopEvent();
#endif
  glutSetWindow(window_id);
  return interrupted;
}

// Screenshot handler
void saveImageAndExit(int dummy)
{
  // set pointmode to 1 to enforce rendering of all visible points
  pointmode = 1;
  if (screenshot_filename.empty()) {
    saveImage(0);
  } else {
    saveImageAt(screenshot_filename);
  }
  exit(0);
}

// GLUT Function wrappers:

void keyPressedIm(unsigned char key, int x, int y)
{
  // We need to update unless we are in an animation
  if (haveToUpdate != 3) haveToUpdate = 1;
  interruptDrawing();
  ImGuiIO &io = ImGui::GetIO();
  ImGui_ImplGLUT_KeyboardFunc(key, x, y);
  if (!io.WantCaptureKeyboard) {
    callbacks::glut::keyPressed(key, x, y);
  }
}

void keyPressedUpIm(unsigned char key, int x, int y)
{
  if (haveToUpdate != 3) haveToUpdate = 1;
  interruptDrawing();
  ImGui_ImplGLUT_KeyboardUpFunc(key, x, y);
  callbacks::glut::keyReleased(key, x, y);
}

/* Wrapper for mouse button handling. Uses both the internal callback and
 * imgui*/
void mouseButtonIm(int button, int state, int x, int y)
{
  if (haveToUpdate != 3) haveToUpdate = 1;
  interruptDrawing();
  ImGui_ImplGLUT_MouseFunc(button, state, x, y);
  callbacks::glut::mouseButton(button, state, x, y);
}

/* Wrapper for reshape handling. Uses both the internal callback and imgui*/
void reshapeIm(int width, int height)
{
  ImGui_ImplGLUT_ReshapeFunc(width, height);
  callbacks::glut::reshape(width, height);
}

/* Wrapper for mouse movement. Uses both the internal callback and imgui*/
void mouseMoveIm(int x, int y)
{
  if (haveToUpdate != 3) haveToUpdate = 1;
  interruptDrawing();
  ImGui_ImplGLUT_MotionFunc(x, y);
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse)
    callbacks::glut::mouseMove(x, y);
}

// Static ImGui state variables
static vgm::Quat qRotVgm = vgm::Quat(1.f, 0.f, 0.f, 0.f);
static double t[3] = {0, 0, 0}, mat[16], rPT[3];
static int colorTypeVal = 0;
static int modal_renderings = 0;
static int colorMapVal = 0;
static bool cam_mouse_nav_bool = true;
static bool always_all_pts = false;
static bool always_reduce_pts = false;
static bool select_unselect_bool = true;
static bool select_voxels_bool = false;
static bool show_points_bool = true;
static bool show_objects_bool = true;
static bool show_cameras_bool = true;
static bool show_path_bool = true;
static bool show_poses_bool = true;
static bool pathshift_bool;
static int signed_int_cam_choice;
static bool save_anim_bool, interpolate_bool;
static bool coloranimbool;

/**
 * Renders ImGui.
 * Includes calls to ImGui::NewFrame() and ImGui::Render().
 * After this function, you have to manually collect and render, e.g.:
 * 	ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
 */
void renderImGuiWindows()
{
  // ImGUI Renderings First:
  ImGui_ImplOpenGL2_NewFrame();
  ImGui_ImplGLUT_NewFrame();
  ImGui::NewFrame();
  ImGuiIO &io = ImGui::GetIO();
  io.ConfigWindowsMoveFromTitleBarOnly = true;

  // 1. Selection Window
  {
    ImGui::SetNextWindowPos(ImVec2(START_WIDTH_IMGUI * 0.83, START_HEIGHT_IMGUI * 0.01),
          ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(START_WIDTH_IMGUI * 0.165, START_HEIGHT_IMGUI * 0.80),
           ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowCollapsed(true, ImGuiCond_FirstUseEver);
    ImGui::Begin("Selection");

    if (ImGui::TreeNode("Draw")) {
      // Drawing Controlls
      // Checkboxes
      ImGui::Checkbox("Draw Points", &show_points_bool);
      ImGui::Checkbox("Draw Objects", &show_objects_bool);
      ImGui::Checkbox("Draw Camera", &show_cameras_bool);
      ImGui::Checkbox("Draw Path", &show_path_bool);
      ImGui::Checkbox("Draw Poses", &show_poses_bool);
      ImGui::Separator();
      // Point size
      ImGui::Text("Point size:");
      ImGui::SliderFloat("##pixels", &pointsize, 0.1f, 10.0f, "%.1f", ImGuiSliderFlags_Logarithmic);
      ImGui::Separator();

      ImGui::TreePop();
    }
    ImGui::Separator();

    // Fog Controlls
    if (ImGui::TreeNode("Fog")) {
      ImGui::RadioButton("Disable", &show_fog, 0);
      ImGui::RadioButton("Exp", &show_fog, 1);
      ImGui::RadioButton("Exp2", &show_fog, 2);
      ImGui::RadioButton("Lin", &show_fog, 3);
      ImGui::RadioButton("Inv Exp", &show_fog, 4);
      ImGui::RadioButton("Inv Exp2 ", &show_fog, 5);
      ImGui::RadioButton("Inv Lin", &show_fog, 6);
      ImGui::Separator();
      ImGui::TreePop();
    }
    ImGui::Separator();
    ImGui::Text("Fog Density:");
    ImGui::SliderFloat("##fog", &fogDensity, 0.0f, 1.0f, "%.5f", ImGuiSliderFlags_Logarithmic);
    ImGui::Separator();

    // Color Controlls
    if (ImGui::TreeNode("Color")) {
      ImGui::Checkbox("Invert", &invert);

      // Color Values
      ImGui::Separator();
      colorTypeVal = listboxColorVal;
      if (ImGui::TreeNode("Values")) {
        ImGui::RadioButton("Height", &listboxColorVal, 0);
        if (!(pointtype.hasReflectance()))
          ImGui::RadioButton("Reflectance", false);
        else
          ImGui::RadioButton("Reflectance", &listboxColorVal, 1);
        if (!(pointtype.hasTemperature()))
          ImGui::RadioButton("Temperature", false);
        else
          ImGui::RadioButton("Temperature", &listboxColorVal, 2);
        if (!(pointtype.hasAmplitude()))
          ImGui::RadioButton("Amplitude", false);
        else
          ImGui::RadioButton("Amplitude", &listboxColorVal, 3);
        if (!(pointtype.hasDeviation()))
          ImGui::RadioButton("Deviation", false);
        else
          ImGui::RadioButton("Deviation", &listboxColorVal, 4);
        if (!(pointtype.hasType()))
          ImGui::RadioButton("Type", false);
        else
          ImGui::RadioButton("Type", &listboxColorVal, 5);
        if (colorTypeVal != listboxColorVal) {
          colorTypeVal = listboxColorVal;
          mapColorToValueIm();
        }
        ImGui::TreePop();
      }
      // Color Map
      ImGui::Separator();
      colorMapVal = listboxColorMapVal;
      if (ImGui::TreeNode("Map")) {
        ImGui::RadioButton("Solid", &listboxColorMapVal, 0);
        ImGui::RadioButton("Grey", &listboxColorMapVal, 1);
        ImGui::RadioButton("HSV", &listboxColorMapVal, 2);
        ImGui::RadioButton("Jet", &listboxColorMapVal, 3);
        ImGui::RadioButton("Hot", &listboxColorMapVal, 4);
        ImGui::RadioButton("Rand", &listboxColorMapVal, 5);
        ImGui::RadioButton("SHSV", &listboxColorMapVal, 6);
        ImGui::RadioButton("TEMP", &listboxColorMapVal, 7);
        if (colorMapVal != listboxColorMapVal) {
          colorMapVal = listboxColorMapVal;
          changeColorMap(0); // dummy 0
        }
        ImGui::TreePop();
      }
      // Color Type
      ImGui::Separator();
      if (ImGui::TreeNode("Type")) {
        ImGui::RadioButton("Use Color Map", &colorScanVal, 0);
        ImGui::RadioButton("By Scan ID", &colorScanVal, 1);
        if (!(pointtype.hasColor()))
          ImGui::RadioButton("Point Color", false);
        else
          ImGui::RadioButton("Point Color", &colorScanVal, 2);
        setScansColored(0); // dummy 0
        ImGui::TreePop();
      }
      ImGui::Separator();
      ImGui::Text("Color values:");
      if (ImGui::Button("Reset Min/Max")) {
        resetMinMax(0);
      }
      ImGui::InputFloat("Min", &mincolor_value, 10, 100, "%.1f", ImGuiInputTextFlags_CharsDecimal);
      ImGui::InputFloat("Max", &maxcolor_value, 10, 100, "%.1f", ImGuiInputTextFlags_CharsDecimal);
      minmaxChanged(0); // dummy0
      ImGui::TreePop();
    }

    // Anim panel
    ImGui::Checkbox("Force keep color", &coloranimbool);
    coloranim = coloranimbool;
    ImGui::Separator();
    ImGui::Text("Anim delay:");
    ImGui::InputInt("##Anim", &anim_delay, 1, 10);
    anim_delay = std::max(0, anim_delay);
    if (ImGui::Button("Animate")) {
      startAnimation(0);
      pointmode = 1; // hack to prevent modal dialog to be rendered this cycle
      // the correct pointmode gets set next cycle in idleIm()
    }
    ImGui::Separator();

    // Camera path
    if (ImGui::TreeNode("Camera Path")) {
      ImGui::Text("Path file:");
      ImGui::InputTextWithHint("##File1", "path.dat", path_file_name, 2048);
      ImGui::Text("Interpolation factor");
      ImGui::SliderFloat("##Interpolationfactor", &path_interp_factor, 0.001f, 10.0f, "%.3f",
             ImGuiSliderFlags_Logarithmic); // Logarithmic slider feels good
      if (ImGui::Button("Save Path")) {
        savePath(0);
      }
      ImGui::SameLine();
      if (ImGui::Button("Load Path")) {
        loadPath(0);
      }
      if (ImGui::Button("Load Robot Path")) {
        drawRobotPath(0);
      }
      ImGui::Separator();
      ImGui::Checkbox("Save Animation", &save_anim_bool);
      save_animation = save_anim_bool;
      ImGui::Checkbox("Interpolate by Distance", &interpolate_bool);
      // Check if updated
      inter_by_dist = interpolate_bool;
      if (ImGui::Button("Animate Path")) {
        updateCamera();
        pathAnimate(0);
        pointmode = 1; // hack to prevent modal dialog to be rendered this cycle
        // the correct pointmode gets set next cycle in idleIm()
      }
      if (ImGui::Button("Animate Path and Matching")) {
        pathMatchingAnimate(0);
        pointmode = 1; // hack to prevent modal dialog to be rendered this cycle
        // the correct pointmode gets set next cycle in idleIm()
      }
      ImGui::TreePop();
    }
    ImGui::Separator();

    // Screenshots and position config
    if (ImGui::TreeNode("Position & Screenshots")) {
      ImGui::Text("Pose file:");
      ImGui::InputTextWithHint("##PoseFile", "pose.dat", pose_file_name, 2048);
      if (ImGui::Button("Save"))
        savePose(0);
      ImGui::SameLine();
      if (ImGui::Button("Load"))
        loadPose(0);
      ImGui::Text("Scaling Factor");
      ImGui::SliderInt("##Factor", &factor, 1, 10);
      if (ImGui::Button("Save Image")) {
        pointmode = 1;
        saveImage(0); // screenshot before imgui frame render x)
      }
      ImGui::TreePop();
    }
    ImGui::Separator();

    // Selection Panel
    if (ImGui::TreeNode("Selection")) {
      ImGui::Text("Selection file:");
      ImGui::InputTextWithHint("##File2", "selected.3d", selection_file_name, 2048);
      if (ImGui::Button("Save"))
        saveSelection(0);
      ImGui::SameLine();
      if (ImGui::Button("Load"))
        loadSelection(0);
      ImGui::SameLine();
      if (ImGui::Button("Clear"))
        clearSelection(0);
      ImGui::Checkbox("Select/Unselect", &select_unselect_bool);
      selectOrunselect = select_unselect_bool;
      ImGui::Checkbox("Select Voxels", &select_voxels_bool);
      select_voxels = select_voxels_bool;
      ImGui::SliderInt("Depth:", &selection_depth, 1, 100, "%d", ImGuiSliderFlags_Logarithmic);
      ImGui::SliderInt("Brushsize:", &brush_size, 0, 100, "%d", ImGuiSliderFlags_Logarithmic);
      ImGui::TreePop();
    }
    ImGui::Separator();

    // Advanced panel
    if (advanced_controls && ImGui::TreeNode("Advanced")) {

      // Scan range selection
      ImGui::Separator();
      ImGui::Text("Scan range selection:");
      ImGui::SliderInt("Id0", &startRangeScanIdx, startScanIdx, endRangeScanIdx);
      ImGui::SliderInt("Idx", &endRangeScanIdx, startRangeScanIdx, endScanIdx);
      if (ImGui::Button("Step Up"))
        stepScansUp(0);
      ImGui::SameLine();
      if (ImGui::Button("Step Down"))
        stepScansDown(0);
      if (ImGui::Button("Reload Frames"))
        reloadFrames();

      ImGui::SliderInt("Frame #:", &current_frame, 0, MetaMatrix[0].size() - 1);
      ImGui::SliderFloat("FPS. Lim:", &idealfps, 0.0f, 240.0f, "%.1f");
      ImGui::SliderFloat("Farplane:", &maxfardistance, 1.0f, 1000000.0f, "%.1f",
             ImGuiSliderFlags_Logarithmic);
      ImGui::SliderFloat("Nearplane:", &neardistance, 0.01f, 1000000.0f, "%.01f",
             ImGuiSliderFlags_Logarithmic);
      if (ImGui::Button("Cylce LOD"))
        cycleLOD();
      ImGui::SliderFloat("LOD speed:", &adaption_rate, 0.0f, 3.0f, "%.1f");
      ImGui::Checkbox("Shift Path for 3D", &pathshift_bool);
      path3D = pathshift_bool;
      ImGui::SliderFloat("3D Shift:", &shifted, 0.0f, 50.0f);
      ImGui::TreePop();

      // idk why this has not been there earlier
    } else if (!advanced_controls) {
      ImGui::Text("--advanced missing");
      ImGui::SameLine();
      if (ImGui::Button("Enable"))
        advanced_controls = true;
      ImGui::Separator();
    }

    ImGui::End(); // End of "3D-Viewer Selection" Panel
  }

  // 2. Similar to legacy shows controlls panel
  {
    ImGui::SetNextWindowPos(ImVec2(START_WIDTH_IMGUI * 0.01, START_HEIGHT_IMGUI * 0.01),
          ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(START_WIDTH_IMGUI * 0.655, START_HEIGHT_IMGUI * 0.20),
           ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowCollapsed(true, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Controls")) {

      bool table_exists = ImGui::BeginTable(
          "controls_table", 5, ImGuiTableFlags_ScrollX | ImGuiTableFlags_SizingStretchProp);

      // Column 1
      ImGui::TableNextColumn();
      ImGui::Text("Mode");
      if (ImGui::Button("Top view"))
        topView();
      if (ImGui::Button("Rotate view"))
        rotateView();
      if (ImGui::Button("Reset view"))
        resetView(0);

      // Column 2
      ImGui::TableNextColumn();
      ImGui::Text("Zoom settings");
      if (showViewMode == 0) {
        ImGui::SliderFloat("Field of View", &cangle, 1.0, 180.0, "%.1f");
        ImGui::Text("Currently in: Normal view.");
      }
      if (showViewMode == 1) {
        ImGui::SliderFloat("Parallel Zoom", &pzoom, 0.1, 100000.0, "%.1f",
               ImGuiSliderFlags_Logarithmic);
        ImGui::Text("Currently in: Ortographic view.");
      }
      if (showViewMode == 2) {
        ImGui::SliderFloat("Rotate Zoom", &rzoom, 0.001, 100000.0, "%.3f",
               ImGuiSliderFlags_Logarithmic);
        ImGui::Text("Currently in: Rotate view.");
      }
      if (ImGui::TreeNode("Keyboard bindings")) {
        ImGui::Text("w/a/s/d: fwd/back left/right");
        ImGui::Text("q/e: roll\t\tc/y: up/down");
        ImGui::Text("Hold Shift/Ctrl: faster/slower");
        ImGui::TreePop();
      }
      // Column 3
      // Navball replacement
      ImGui::TableNextColumn();
      ImGui::Text("Navball");
      if (ImGui::gizmo3D("##gizmo1", qRotVgm)) {
        if (showViewMode != 2) {
          // Convert gizmos vgm crap to 3dtk
          quat[0] = qRotVgm.w;
          quat[1] = qRotVgm.x;
          quat[2] = qRotVgm.y;
          quat[3] = qRotVgm.z;
          QuatToMatrix4((const double *)quat, t, mat);
          Matrix4ToEuler(mat, rPT);
          // rotate the camera
          mouseRotX = deg(rPT[0]);
          mouseRotY = deg(rPT[1]);
          mouseRotZ = deg(rPT[2]);
        }
      } else {
        rPT[0] = rad(mouseRotX);
        rPT[1] = rad(mouseRotY);
        rPT[2] = rad(mouseRotZ);
        EulerToMatrix4(t, rPT, mat);
        // Convert back to vgm, which gets synched next
        // rendering cycle
        if (showViewMode != 2) {
          Matrix4ToQuat(mat, quat);
          qRotVgm.w = quat[0];
          qRotVgm.x = quat[1];
          qRotVgm.y = quat[2];
          qRotVgm.z = quat[3];
        }
      }

      // Column 4
      ImGui::TableNextColumn();
      ImGui::Text("Camera");
      ImGui::SliderInt("Choose camera", &signed_int_cam_choice, 0, cams.size());
      cam_choice = signed_int_cam_choice;
      if (ImGui::Button("Add camera"))
        callAddCamera(0);
      if (ImGui::Button("Delete camera"))
        callDeleteCamera(0);

      // Column 5
      ImGui::TableNextColumn();
      ImGui::Text("General behaviour");
      ImGui::Checkbox("MouseNav", &cam_mouse_nav_bool);
      bool always_tmp = always_all_pts, reduce_tmp = always_reduce_pts;
      ImGui::Checkbox("Always all points", &always_all_pts);
      ImGui::Checkbox("Always reduce points", &always_reduce_pts);

      // There might only be one checkbox marked at a time:

      // Check if changed always pts
      if (always_tmp != always_all_pts && always_all_pts) {
        always_reduce_pts = false;
        reduce_tmp = false;
      } // Check if changed reduced pts
      else if (reduce_tmp != always_reduce_pts && always_reduce_pts) {
        always_all_pts = false;
        always_tmp = false;
      }

      // Close table. But only if the screen is large enough
      // such that is exists
      if (table_exists)
        ImGui::EndTable();
    }
    ImGui::End(); // End of "3D Viewer - Controlls panel"
  }

  /**
   * Insanity of Show states (How to decide pointmode):
   *
   * pointmode is an int that has following modes:
   *  pointmode == -1 : always reduce points
   *  pointmode ==  0 : accumulate points to screen directly
   *  pointmode ==  1 : always render all points
   *
   * These internal modes must be decided from the ImGui checkboxes
   * always_reduce_pts and always_all_pts. However, when interacted
   * with the scene, pointmode should switch to -1 (reduce all), and
   * must then be restored when no longer interacted.
   **/

  // If Idle and no checkbox marked:
  if (!mousemoving && !keypressed && !always_reduce_pts && !always_all_pts) {
    // Change pointmode to display everything
    if (pointmode != 0) {
      pointmode = 0;
    }

    // If non-idle and always reduce
  } else if (always_reduce_pts) {
    // Change pointmode to always reduce
    if (pointmode != -1) {
      pointmode = -1;
    }

    // If non-idle and always all
  } else if (always_all_pts) {
    // interuptable always_all_pts mode
    checkForInterrupt();
    // Start loading the points, show a modal dialog (render twice because rendering cycles...)
    if (!isInterrupted() && pointmode != 1 && !mousemoving && !keypressed && modal_renderings < 2 &&
        haveToUpdate != 6) {
      ImGui::OpenPopup("Please wait");
      if (ImGui::BeginPopupModal("Please wait")) {
        ImGui::Text("Loading all the points at once...");
        ImGui::Text("Interact to abort.");
        ImGui::EndPopup();
      }
      modal_renderings++;
      haveToUpdate = 1;
      // If we are interrupted, reset mode to display only reduced points internally
    } else if (mousemoving || keypressed || isInterrupted()) {
      if (pointmode != -1) {
        pointmode = -1;
        if (haveToUpdate != 3)
          haveToUpdate = 1;
      }
      if (ImGui::IsPopupOpen("Please wait")) {
        ImGui::CloseCurrentPopup();
      }
      modal_renderings = 0;

      // If we did not get interrupted and have rendered the modal dialog, start rendering all points
    } else if (modal_renderings >= 2 && pointmode != 1) {
      // Change pointmode to always all
      pointmode = 1;
      if (ImGui::IsPopupOpen("Please wait")) {
        ImGui::CloseCurrentPopup();
      }
      modal_renderings = 0;
      glutPostRedisplay();
    }
    // If no idle and no checkbox
  } else if (!always_reduce_pts && !always_all_pts && (mousemoving || keypressed)) {
    // Change to always reduce
    if (pointmode != -1) {
      pointmode = -1;
      glutPostRedisplay();
    }
  }

  // Set show states to ImGui states
  cameraNavMouseMode = cam_mouse_nav_bool;
  show_points = show_points_bool;
  show_objects = show_objects_bool;
  show_cameras = show_cameras_bool;
  show_path = show_path_bool;
  show_poses = show_poses_bool;

  // Render IMGUI stuff into ImGuis internal buffer
  ImGui::Render();
}

/**
 * Displays all the points.
 * This function is similar to DrawPoints in show_gl.cc.
 * However, this function will also draw ImGui directly to GL_FRONT if interruptable is true.
 *
 * @param mode: Usually GL_RENDERER
 * @param interruptable: If true, render directly to the screen and not to the back buffer GL_BACK.
 */
void DrawPointsIm(GLenum mode, bool interruptable)
{
  long time = GetCurrentTimeInMilliSec();
  double min = 0.000000001;
  double max = 1.0;
  LevelOfDetail *= 1.0 + adaption_rate * (lastfps - idealfps) / idealfps;
  if (LevelOfDetail > max)
    LevelOfDetail = max;
  else if (LevelOfDetail < min)
    LevelOfDetail = min;

  // In case of animation
  if (frameNr != 0) {
    if (coloranim == 0) {
      cm->setMode(ScanColorManager::MODE_ANIMATION);
    }

    for (int iterator = (int)octpts.size() - 1; iterator >= 0; iterator--) {

      // ignore scans that don't have any frames associated
      // with them
      if ((unsigned int)iterator >= MetaMatrix.size())
        continue;

      // also ignore scans outside the selected range - if in
      // advanced mode
      if (advanced_controls) {
        // pay attention to offset (startScanIdx)
        if (iterator < startRangeScanIdx - startScanIdx)
          continue;
        if (iterator > endRangeScanIdx - startScanIdx)
          continue;
      }

      // set usable frame
      double *frame;
      Scan::AlgoType type;
      if ((unsigned int)frameNr >= MetaMatrix[iterator].size()) {
        // use last possible frame
        frame = MetaMatrix[iterator].back();
        type = MetaAlgoType[iterator].back();
      } else {
        frame = MetaMatrix[iterator][frameNr];
        type = MetaAlgoType[iterator][frameNr];
      }
      if (type == Scan::INVALID)
        continue;
      cm->selectColors(type);
      glPushMatrix();
      glMultMatrixd(frame);
      glPointSize(pointsize);
      ExtractFrustum(pointsize);
      cm->selectColors(type);
      checkForInterrupt();
      if (pointmode == 1 && !mousemoving && !keypressed && !isInterrupted()) {
        octpts[iterator]->display();
      } else {
        octpts[iterator]->displayLOD(LevelOfDetail);
      }
      glPopMatrix();
    }

    setScansColored(0);

  } else {

    if (mode == GL_SELECT) {
      // select points mode
      // ------------------
      GLuint name = 0;
      for (int iterator = (int)octpts.size() - 1; iterator >= 0; iterator--) {
        glPushMatrix();
        glMultMatrixd(MetaMatrix[iterator].back());

        glColor4f(1.0, 0.0, 0.0, 1.0);
        glPointSize(pointsize + 2.0);
        for (std::set<sfloat *>::iterator it = selected_points[iterator].begin();
             it != selected_points[iterator].end(); it++) {
          glLoadName(name++);
          glBegin(GL_POINTS);
          glVertex3d((*it)[0], (*it)[1], (*it)[2]);
          glEnd();
        }
        glPointSize(pointsize);

        glFlush();
        glPopMatrix();
      }

    } else {

      // draw point is normal mode
      // -------------------------

      if (interruptable) {
        glDrawBuffer(GL_FRONT);
      }
      glPointSize(pointsize);

      std::vector<int> sequence;
      calcPointSequence(sequence, current_frame);
      for (unsigned int i = 0; i < sequence.size(); i++) {
        int iterator = sequence[i];
        // ignore scans that don't have any frames
        // associated with them
        if ((unsigned int)iterator >= MetaMatrix.size())
          continue;

        // also ignore scans outside the selected range
        // - if in advanced mode
        if (advanced_controls) {
          // pay attention to offset
          // (startScanIdx)
          if (iterator < startRangeScanIdx - startScanIdx)
            continue;
          if (iterator > endRangeScanIdx - startScanIdx)
            continue;
        }

        // set usable frame
        double *frame;
        Scan::AlgoType type;
        if ((unsigned int)current_frame >= MetaMatrix[iterator].size()) {
          // use last possible frame
          frame = MetaMatrix[iterator].back();
          type = MetaAlgoType[iterator].back();
        } else {
          frame = MetaMatrix[iterator][current_frame];
          type = MetaAlgoType[iterator][current_frame];
        }
        if (type == Scan::INVALID)
          continue;
        glPushMatrix();
        if (invert)
          // default: white points on black
          // background
          glColor4d(1.0, 1.0, 1.0, 0.0);
        else
          // black points on white background
          glColor4d(0.0, 0.0, 0.0, 0.0);

        // glMultMatrixd(MetaMatrix[iterator].back());
        if (current_frame != (int)MetaMatrix.back().size() - 1) {
          if (coloranim == 0) {
            cm->setMode(ScanColorManager::MODE_ANIMATION);
            cm->selectColors(type);
          } else {
            setScansColored(0);
          }
        }
        glMultMatrixd(frame);

        ExtractFrustum(pointsize);
        checkForInterrupt();
        if (pointmode == 1 && !mousemoving && !keypressed && !isInterrupted()) {
          octpts[iterator]->display();
        } else if (interruptable) {
          checkForInterrupt();
          // ATTENTION: We sneak ImGui here when
          // force-drawing the points ImGui
          // Windows should never be
          // color-inverted:
          if (!invert) {
            glDisable(GL_COLOR_LOGIC_OP);
          }
          // Sneak imgui
          renderImGuiWindows();
          ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
          // Reset color to before state
          if (invert) {
            glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], 1.0);
          } else {
            glEnable(GL_COLOR_LOGIC_OP);
            glLogicOp(GL_COPY_INVERTED);
            glClearColor(1 - bgcolor[0], 1 - bgcolor[1], 1 - bgcolor[2], 1.0);
          }
          // Force draw everything
          glFlush();
          glFinish();
          if (isInterrupted()) {
            glPopMatrix();
            return;
          }
          octpts[iterator]->display();
        } else {
          octpts[iterator]->displayLOD(LevelOfDetail);
        }
        if (!selected_points[iterator].empty()) {
          glColor4f(1.0, 0.0, 0.0, 1.0);
          glPointSize(pointsize + 2.0);
          glBegin(GL_POINTS);
          for (std::set<sfloat *>::iterator it = selected_points[iterator].begin();
               it != selected_points[iterator].end(); it++) {
            glVertex3d((*it)[0], (*it)[1], (*it)[2]);
          }
          glEnd();
          glPointSize(pointsize);
        }

        glPopMatrix();
      }
    }
  }

  checkForInterrupt();
  if (pointmode == 1 && !mousemoving && !keypressed && !isInterrupted()) {
    fullydisplayed = true;
  } else {
    unsigned long td = (GetCurrentTimeInMilliSec() - time);
    if (td > 0)
      lastfps = 1000.0 / td;
    else
      lastfps = 1000.0;
    fullydisplayed = false;
  }
  if (interruptable)
    fullydisplayed = true;
}

/**
 * Displays the scene (points, poses, path, cameras, etc.)
 * This function is similar to DisplayItFunc in show_gl.cc.
 * However, this function uses DrawPointsIm instead of DrawPoints
 * to support ImGui rendering while rendering to GL_FRONT directly.
 *
 * @param mode: Usually GL_RENDERER
 * @param interruptable: If true, render points directly to the screen
 * and not to the back buffer GL_BACK.
 */
void DisplayItFuncIm(GLenum mode, bool interruptable)
{
  // set the clear color buffer in case of
  // both invert and non invert mode
  if (invert) {
    glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], 1.0);
  } else {
    glEnable(GL_COLOR_LOGIC_OP);
    glLogicOp(GL_COPY_INVERTED);
    glClearColor(1 - bgcolor[0], 1 - bgcolor[1], 1 - bgcolor[2], 1.0);
  }

  // clear the color and depth buffer bit
  if (!interruptable) { // single buffer mode, we need the depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  }

  // set the polygon mode XXX please clarify why we need this and what GL_FRONT_AND_BACK would do
  glPolygonMode(GL_FRONT /*_AND_BACK*/, GL_LINE);

  glPushMatrix();

  setup_camera();

  DrawScala();

  setup_fog();

  if (fardistance > maxfardistance)
    fardistance = maxfardistance;
  if (fabs(oldfardistance - fardistance) > 0.00001 || fabs(oldneardistance - neardistance) > 0.00001) {
    oldfardistance = fardistance;
    oldneardistance = neardistance;
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    reshapeIm(viewport[2], viewport[3]);
  }

  if (show_path == 1) {
    double *pose;
    glColor4d(1.0, 0.0, 0.0, 1.0);
    glLineWidth(5);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; (unsigned int)i < MetaMatrix.size(); i++) {

      // also ignore scans outside the selected range - if in
      // advanced mode
      if (advanced_controls) {
        // pay attention to offset (startScanIdx)
        if (i < startRangeScanIdx - startScanIdx)
          continue;
        if (i > endRangeScanIdx - startScanIdx)
          continue;
      }
      // set usable type
      Scan::AlgoType type;
      if ((unsigned int)frameNr >= MetaMatrix[i].size()) {
        type = MetaAlgoType[i].back();
      } else {
        type = MetaAlgoType[i][frameNr];
      }
      if (frameNr >= 1 && frameNr < (int)MetaMatrix[i].size()) {
        if (type == Scan::INVALID)
          continue;
        // avoid incomplete frames in a scan
        if ((unsigned int)frameNr >= MetaMatrix[i].size())
          pose = MetaMatrix[i].back();
        else
          pose = MetaMatrix[i][frameNr];
      } else {
        //  avoid incomplete frames in a scan
        if ((unsigned int)current_frame >= MetaMatrix[i].size())
          pose = MetaMatrix[i].back();
        else
          pose = MetaMatrix[i][current_frame];
      }
      if (showViewMode == 1) {
        glVertex3f(pose[12], pose[13], pose[14]);
      } else {
        glVertex3f(pose[12], pose[13], pose[14]);
      }
    }
    glEnd();
  }

  // Draw trajectory from file
  if (trajectory.size() > 0) {
    glColor4d(1.0, 0.0, 0.0, 1.0);
    glLineWidth(5);
    glBegin(GL_LINE_STRIP);

    for (size_t i = 0; i < trajectory.size(); i++) {
      glVertex3f(trajectory.at(i)[0], trajectory.at(i)[1], trajectory.at(i)[2]);
    }

    glEnd();
  }

  // if show camera is true then draw cameras.
  if (show_cameras == 1) {
    DrawCameras();
  }

  // if show path is true the draw path.
  if (show_path == 1) {
    DrawPath();
  }
  // if show poses is true then draw coordinate axes.
  if (show_poses == 1) {
    DrawCoordinateSystems();
  }
  // if show objects is true then draw objects.
  if (show_objects == 1) {
    DrawObjects(mode);
  }

  // if show points is true then draw points
  if (show_points == 1)
    DrawPointsIm(mode, interruptable);

  // if show_cylinderBody or show_cylinderPoints true draw cylinder
  if (show_cylinderBody || show_cylinderPoints)
    DrawCylinder();

  if (classLabels)
    DrawTypeLegend();
  if (label)
    DrawUrl();

  glPopMatrix();

  if (!invert) {
    glDisable(GL_COLOR_LOGIC_OP);
  }
}

/**
 * Display loop function for GLUT.
 * This function gets called the next rendering cycle if glutPostRedisplay() has been called.
 * If glutPostRedisplay() has not been called, the idleIm() function is called instead.
 * This function is responsible for swapping the GL_BACK and GL_FRONT buffers, i.e.,
 * it calls glutSwapBuffers(), which idleIm() does not call.
 */
void displayIm()
{
  // Render ImGui windows
  renderImGuiWindows();

  // GLUT / OpenGL2 camera and aspect handling:
  if (((fabs(cangle_old - cangle) > 0.5)) || (fabs(pzoom_old - pzoom) > 0.5)) {
    cangle_old = cangle;
    pzoom_old = pzoom;
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    reshapeIm(viewport[2], viewport[3]);
#ifdef _WIN32
    Sleep(5); // legacy show: Sleep(25)
#else
    usleep(50000); // legacy show: usleep(250000)
#endif
  }

  // show the rendered scene
  glDrawBuffer(buffermode);
  DisplayItFuncIm(GL_RENDER, false);

  // Finally Draw ImGui contents as well
  ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

  // Draw the buffer
  glutSwapBuffers();
  // Update if ImGui wants attention or reduce pts (interactive mode)
  /* TODO : Think about a fourth mode, similar to -1 (reduce all), but only when interacted */
  if (pointmode == -1 || ImGui::GetIO().WantCaptureMouse || ImGui::GetIO().WantCaptureKeyboard) {
    if (pointmode != 0)
      glutPostRedisplay();
  }
}

/**
 * Idle function for GLUT.
 * This function gets called every rendering cycle if there are no pending events or glutPostRedisplay() has not been
 * called. If glutPostRedisplay() has been called or some other event triggers a redraw, the displayIm() function is
 * called instead. This function is not responsible for swapping the GL_BACK and GL_FRONT buffers, i.e., it does not
 * call glutSwapBuffers(), which displayIm() calls.
 */
void idleIm(void)
{
#ifdef _WIN32
  Sleep(1);
#else
  usleep(1000);
#endif


  if (glutGetWindow() != window_id)
    glutSetWindow(window_id);

  // return as nothing has to be updated
  if (haveToUpdate == 0) {
    if (!fullydisplayed && !mousemoving && !keypressed && pointmode == 0) {
      glDrawBuffer(buffermode);
      // Call the display function
      DisplayItFuncIm(GL_RENDER, true);
    }
    return;
  }

  // case: display is invalid - update it
  if (haveToUpdate == 1) {
    update_callback();
    haveToUpdate = 0;
    return;
  }

  // case: camera angle is changed - instead of repeating code call
  // Reshape, since these OpenGL commands are the same
  if (haveToUpdate == 2) {
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    reshapeIm(viewport[2], viewport[3]);
    update_callback();
    haveToUpdate = 0;
    return;
  }

  // case: animation
  if (haveToUpdate == 3) {
    frameNr += 1;
    if (!(MetaMatrix.size() > 1 && frameNr < (int)MetaMatrix[1].size())) {
      frameNr = 0;
      haveToUpdate = 4;
      return;
    }
    update_callback();

    if (save_animation) {
      std::string filename = scan_dir + "animframe" + to_string(frameNr, 5) + ".png";
      std::cout << "write " << filename << std::endl;
      int tmpUpdate = haveToUpdate;
      glWriteImagePNG(filename.c_str(), factor, 0);
      haveToUpdate = tmpUpdate;
    }
  }

#ifdef _WIN32
  // Extra Sleep is NOT necessary!
  // Sleep(300);
  Sleep(anim_delay);
#else
  usleep(anim_delay * 10000);
#endif

  if (haveToUpdate == 4) { // stop animation
    frameNr = 0;	 // delete these lines if you want a 'continue' functionality.
    haveToUpdate = 1;
  }

  // case: scan matching and path animation in lock-step
  if (haveToUpdate == 8) {
    if (path_iterator == 0) {
      oldcamNavMode = cameraNavMouseMode; // remember state of old mousenav
      cameraNavMouseMode = 0;
    }
    frameNr += 1;
    if (!(MetaMatrix.size() > 1 && frameNr < (int)MetaMatrix[1].size())) {
      frameNr = 0;
      haveToUpdate = 4;
      return;
    }
    std::cout << path_iterator << " " << ups.size() << std::endl;
    if ((false && path_iterator < path_vectorX.size()) ||
        (true && path_iterator < ups.size())) { // standard animation case

      // call the path animation function
      // hide both the cameras and the path
      show_cameras = 0;
      show_path = 0;
      show_poses = 0;

      // increase the iteration count
      path_iterator += 1;
      // repaint the screen
      update_callback();

      // save the animation
      if (save_animation) {
        std::string filename = scan_dir + "animframe" + to_string(path_iterator, 5) + ".png";
        std::cout << "written " << filename << " of " << path_vectorX.size() << " files"
            << std::endl;
        glWriteImagePNG(filename.c_str(), factor, 0);
        haveToUpdate = 8;
      }
    } else { // animation has just ended
      cameraNavMouseMode = oldcamNavMode;
      show_cameras = 1;
      show_path = 1;
      show_poses = 1;
      haveToUpdate = 0;
    }
  }

  // case: path animation
  if (haveToUpdate == 6) {

    if (path_iterator == 0) {
      oldcamNavMode = cameraNavMouseMode; // remember state of
      // old mousenav
      cameraNavMouseMode = 0;
    }

    // check if the user wants to animate both
    // scan matching and the path at the same
    // time
    if (path_iterator < path_vectorX.size()) { // standard animation case

      // call the path animation function
      // hide both the cameras and the path
      show_cameras = 0;
      show_path = 0;
      // increase the iteration count

      path_iterator += 1;
      // repaint the screen
      update_callback();

      // save the animation
      if (save_animation) {
        std::string filename = scan_dir + "animframe" + to_string(path_iterator, 5) + ".png";

        std::cout << "written " << filename << " of " << path_vectorX.size() << " files"
            << std::endl;
        glWriteImagePNG(filename.c_str(), factor, 0);
        haveToUpdate = 6;
      }
    } else { // animation has just ended
      cameraNavMouseMode = oldcamNavMode;
      show_cameras = 1;
      show_path = 1;
      haveToUpdate = 0;
    }
  }
}

/**
 * This function installs the callback functions for GLUT (including ImGui wrappers)
 * and determines initial window position, width, and heigth.
 */
void initScreenWindowIm()
{
  // init display
  glutInitDisplayMode(GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE);
  SCREEN_WIDTH = glutGet(GLUT_SCREEN_WIDTH); // This works well, when
  // there is only one monitor.
  SCREEN_HEIGHT = glutGet(GLUT_SCREEN_HEIGHT); // For multimonitor this
  // needs more treatment...
  // Stick to original 3dtk
  START_WIDTH_IMGUI = 1280; // Factor of 2 will give QHD, Andreas ;)
  START_HEIGHT_IMGUI = 720;
  // Check if screen is large enough for HD...
  if (START_WIDTH_IMGUI > SCREEN_WIDTH || START_HEIGHT_IMGUI > SCREEN_HEIGHT) {
    // Use original 3dtk
    START_WIDTH_IMGUI = START_WIDTH;
    START_HEIGHT_IMGUI = START_HEIGHT;
  }
  BOOST_ASSERT_MSG(!(START_WIDTH_IMGUI > SCREEN_WIDTH || START_HEIGHT_IMGUI > SCREEN_HEIGHT),
       "Go get a wider screen.");
  START_X_IMGUI = 0; //(SCREEN_WIDTH - START_WIDTH_IMGUI) / 2;
  START_Y_IMGUI = 0; //(SCREEN_HEIGHT - START_HEIGHT_IMGUI) / 2;

  // define the window position and size
  glutInitWindowPosition(START_X_IMGUI, START_Y_IMGUI);
  glutInitWindowSize(START_WIDTH_IMGUI, START_HEIGHT_IMGUI);

  // create window
  window_id = glutCreateWindow("3D_Viewer");
  // Set the loop function for GLUT
  glutDisplayFunc(displayIm);

  // Initialize ImGui, get context and backend
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui::StyleColorsDark();

  // Disable imgui.ini file which messes with initial positions sometimes
  ImGui::GetIO().IniFilename = nullptr;
  ImGui_ImplGLUT_Init();
  ImGui_ImplOpenGL2_Init();

  // Install GLUT functions
  // ImGui wraps the original 3dtk GLUT handlers:
  glutReshapeFunc(reshapeIm);
  glutMouseFunc(mouseButtonIm);
  glutMotionFunc(mouseMoveIm);
#ifdef __FREEGLUT_EXT_H__
  glutMouseWheelFunc(ImGui_ImplGLUT_MouseWheelFunc);
#endif
  // Keyboard and Idle are not handled by 3dtk.
  // Original 3dtk installed GLUT handler functions:
  glutKeyboardFunc(keyPressedIm);
  glutKeyboardUpFunc(keyPressedUpIm);
  // glutSpecialFunc(ImGui_ImplGLUT_SpecialFunc);
  glutIdleFunc(idleIm);
#ifdef __APPLE__
  glutWMCloseFunc(quit); // apple works different
#else
  glutCloseFunc(callbacks::glut::quit);
#endif

  // Setup view and 3DTK logo
  update_view_rotate(0);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  load_url_texture();
}

/**
 * Main function.
 * Reads the scan (scan000.3d, ...) and frames files (scan000.frames, ...)
 * from the data directory.
 * The frames are used for animation of the matching process.
 */
int main(int argc, char **argv)
{
  // Signal handling and GLUT init
  setSignalHandling();
  glutInit(&argc, argv);

  // Program options
  dataset_settings dss;
  window_settings ws;
  display_settings ds;
  try {
    parse_args(argc, argv, dss, ws, ds);
  } catch (std::exception &e) {
    std::cerr << "Error while parsing settings: " << e.what() << std::endl;
    exit(1);
  }

  // Connect backend and init show window
  initScreenWindowIm();
  initShow(dss, ws, ds);
  // Start in accumulation mode and disable the bottom left label
  pointmode = 0;
  label = false;
  if (takescreenshot) {
    glutTimerFunc(0, &saveImageAndExit, 0);
  }
#ifndef __APPLE__
  std::thread t1(spacenavHandlerIm);
#endif

  // Starting the Glut Main Loop
  glutMainLoop();

  // Cleanup
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGLUT_Shutdown();
  ImGui::DestroyContext();
}
