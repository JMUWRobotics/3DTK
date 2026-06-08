New Click Tool:

New Click Tool is a subproject of the 3DTK project that allows users to select and export points in a single 2D image or corresponding points between two 2D images either loaded directly as images or as 3D scans which are internally converted to 2D images. It provides a zoom function to select points with an accuracy of up to 0.1 pixels. 


Structure:
  
    3DTK/
    |--src/slam6d/fbr/new_click_tool/
    | |--docs/
    | | |--Benutzerhandbuch.pdf
    | | |--Entwicklerhandbuch.pdf
    | | |--Pflichtenheft.pdf
    | |--MVP /
    | | |--assets/
    | | | |--mein_Bild.jpg
    | | | |--...
    | | |--CMakeLists.txt
    | | |--new_click_tool_app.cpp
    | | |--new_click_tool_app_two.cpp
    | | |--new_click_tool_main.cpp
    | | |--...
    |--include/slam6d/fbr
    | |--...
    | |--new_click_tool_app.h
    | |--...
    |--bin/
    | |--...
    | |--new_click_tool
    | |--...


Dependencies:

- C++ standard library
- stb_image
- Dear ImGui
- GLFW
- OpenGL


Build:

The tool is built with the main project (see README.md of the 3DTK project) or as stand-alone tool.
Stand-alone build:
Do the steps of the main project build until the CMake step. After the cmake command, run make new_click_tool instead of make.

Run:

The tool can be started via the binary in the 3DTK/bin/ folder or via command line.

    Some command line parameters: 
    ./new_click_tool                                                   open tool without images
    ./new_click_tool -help                                             shows command line help 
    ./new_click_tool -im <path to image>                               open with one image loaded 
    ./new_click_tool -scan <path to scan>                              open with one scan loaded 
    ./new_click_tool -im <path to image 1> -im <path to image 2>       open with two images loaded 
    ./new_click_tool -im <path to image> -out <path to directory>      open with one image loaded and manually defined folder to save converted panoramas and coordinates

  


