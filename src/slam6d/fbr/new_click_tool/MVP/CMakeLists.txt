cmake_minimum_required(VERSION 3.14)

#ImGUI
# from fbr-CMakesList changed to OpenGL3 and glfw
file(GLOB IMGUI_SRCS_CLICK 
"${PROJECT_SOURCE_DIR}/3rdparty/imgui/*.cpp" 
"${PROJECT_SOURCE_DIR}/3rdparty/imgui/backends/imgui_impl_opengl3.cpp"
"${PROJECT_SOURCE_DIR}/3rdparty/imgui/backends/imgui_impl_glfw.cpp")

add_library(imgui STATIC
    ${IMGUI_SRCS_CLICK}
)

target_include_directories(imgui PUBLIC
    ${IMGUI_INCLUDE_DIRS}
)

add_library(glad STATIC
    ${PROJECT_SOURCE_DIR}/3rdparty/glad/src/glad.c
)
# glad
target_include_directories(glad PUBLIC
    ${PROJECT_SOURCE_DIR}/3rdparty/glad/include
)

add_executable(new_click_tool
    new_click_tool_main.cpp
    new_click_tool_app.cpp
    new_click_tool_app_two.cpp
)

target_include_directories(new_click_tool PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/3rdparty/stb
    ${PROJECT_SOURCE_DIR}/3rdparty/imgui
 #   ${PROJECT_SOURCE_DIR}/3rdparty/imgui/backends
    ${PROJECT_SOURCE_DIR}/3rdparty/glad/include
    ${OpenCV_INCLUDE_DIRS}
)

target_link_libraries(new_click_tool PRIVATE
    glfw
    imgui
    glad
    #for scan-conversion
    scan_to_panorama_lib
    ${OpenCV_LIBS}
    fbr
    scan
    show
    OpenMP::OpenMP_CXX
    
)

set_target_properties(new_click_tool PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)
