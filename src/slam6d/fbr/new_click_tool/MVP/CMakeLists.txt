cmake_minimum_required(VERSION 3.15)
project(ImGuiTool)

set(CMAKE_CXX_STANDARD 17)


#find_package(imgui CONFIG REQUIRED)
#find_package(glfw3 CONFIG REQUIRED)
#find_package(glad CONFIG REQUIRED)
#find_package(glfw3 REQUIRED)


find_package(OpenGL REQUIRED) 

#add_executable(new_click_tool main.cpp App.cpp)
#add_executable(new_click_tool main.cpp App.cpp ${IMGUI_SRCS}) ##
add_executable(new_click_tool main.cpp App.cpp

${CMAKE_SOURCE_DIR}/3rdparty/imgui/imgui.cpp
${CMAKE_SOURCE_DIR}/3rdparty/imgui/imgui_draw.cpp
${CMAKE_SOURCE_DIR}/3rdparty/imgui/imgui_tables.cpp
${CMAKE_SOURCE_DIR}/3rdparty/imgui/imgui_widgets.cpp
${CMAKE_SOURCE_DIR}/3rdparty/imgui/backends/imgui_impl_glfw.cpp
${CMAKE_SOURCE_DIR}/3rdparty/imgui/backends/imgui_impl_opengl3.cpp)

target_include_directories(new_click_tool PRIVATE ${IMGUI_INCLUDE_DIRS}) ##

target_link_libraries(new_click_tool PRIVATE
   # imgui::imgui
    glfw
    #glad::glad
    OpenGL::GL   )

add_compile_definitions(GL_SILENCE_DEPRECATION)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
