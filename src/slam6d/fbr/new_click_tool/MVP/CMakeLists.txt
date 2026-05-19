cmake_minimum_required(VERSION 3.14)
project(ImGuiKoorTool CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ==========================================
# Lokale Bibliotheken aus 3rdparty/
# ==========================================
set(THIRDPARTY_DIR ${CMAKE_SOURCE_DIR}/../../../../../3rdparty)

# 1. GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(${THIRDPARTY_DIR}/glfw ${CMAKE_CURRENT_BINARY_DIR}/glfw_build)

# 2. ImGui
set(IMGUI_DIR ${THIRDPARTY_DIR}/imgui)
add_library(imgui STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC 
    ${IMGUI_DIR} 
    ${IMGUI_DIR}/backends
)
target_link_libraries(imgui PUBLIC glfw)

# 3. GLAD
set(GLAD_DIR ${THIRDPARTY_DIR}/glad)
add_library(glad STATIC ${GLAD_DIR}/src/glad.c)
target_include_directories(glad PUBLIC ${GLAD_DIR}/include)

# ==========================================
# Haupt-Executable & Linken
# ==========================================
add_executable(ImGuiKoorTool main.cpp App.cpp)

# 4. stb_image (nur Header, daher direkt ans Executable gehangen)
target_include_directories(ImGuiKoorTool PRIVATE ${THIRDPARTY_DIR}/stb)

# Linken der Bibliotheken
target_link_libraries(ImGuiKoorTool PRIVATE imgui glad glfw)

# macOS Spezifische Frameworks
if (APPLE)
    find_library(COCOA_FRAMEWORK Cocoa)
    find_library(OPENGL_FRAMEWORK OpenGL)
    find_library(IOKIT_FRAMEWORK IOKit)
    find_library(COREVIDEO_FRAMEWORK CoreVideo)
    
    target_link_libraries(ImGuiKoorTool PRIVATE
        ${COCOA_FRAMEWORK}
        ${OPENGL_FRAMEWORK}
        ${IOKIT_FRAMEWORK}
        ${COREVIDEO_FRAMEWORK}
    )
endif()
