// #define GL_SILENCE_DEPRECATION

// #include <glad/glad.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <filesystem>
#include <iostream>

#include "slam6d/fbr/new_click_tool_app.h"

void usage(char **argv)
{
	std::cout << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << std::endl;
	std::cout << "-help \t \t \t \t \t \t \t Option list" << std::endl;
    std::cout << std::endl;
	std::cout << " No Arguments \t \t \t \t \t \t Open empty window" << std::endl;
	std::cout << "-im <imagePath> \t \t \t \t \t Open one image, default directory" << std::endl;
	std::cout << "-scan <scanPath> \t \t \t \t \t Open one image based on 3D-scan, default directory" << std::endl;
	std::cout << "-scan <scanPath> -im <imagePath> \t \t \t Open two images, default directory" << std::endl;
	std::cout << "-im <imagePath> ... -out <outDir> \t \t \t Set output Directory" << std::endl;
	std::cout << std::endl;
	std::cout << "Default output directory = directory of first loaded image/scan" << std::endl;
    std::cout << "Converted scans: output directory or directory of scan" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "imFormat:" << std::endl;
	std::cout << "-im \t \t \t 2D-image Format" << std::endl;
	std::cout << "-scan \t \t \t 3D-scan Format, internally converted with scan_to_panorama" << std::endl;
	std::cout << std::endl;
	std::cout << "Example: \t \t -im /pictures/myPicture.png -scan /scans/scan001.3D -out /users/desktop/click_tool \n"
		  << std::endl;
	std::cout << std::endl;
}




int main(int argc, char **argv)
{
	if (!glfwInit())
		return 1;

	const char *glsl_version;
#if defined(__APPLE__)
	glsl_version = "#version 150";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#else
	glsl_version = "#version 130";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
#endif

    // Instantiate app
    	App app;

	GLFWwindow *window = glfwCreateWindow(1280, 720, "New click tool", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return 1;
	}
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);
	/*
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
	return -1;
	}*/

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	// find initial image (commandline or emptywindow)
	std::string startImage = "";
	std::string startImage2 = "";
	std::string outputDir = "";
	std::string startScan = "";
	std::string startScan2 = "";
	std::string scanDir = "";
	std::string scanDir2 = "";
    bool twoImageMode = false;

	if (argc > 1) { // if no picture is loaded start app without initial images
		for (int i = 1; i < argc; i++) {
			std::string argString = argv[i];

			if (argString == "-im") { // 2D-image Format
				if (i + 1 < argc) {
					if (startImage.empty() && startScan.empty()){
						startImage = argv[++i];
                        if(!std::filesystem::exists(startImage)){
                            std::cout << "\033[31m Cannot find image " <<std::filesystem::path(startImage).filename().string() << "\033[0m" << std::endl;
                            usage(argv);
                        return 1;
                        }}
					else if (startImage2.empty() && startScan2.empty()){
                        twoImageMode = true;
                        startImage2 = argv[++i];
                        if(!std::filesystem::exists(startImage2)){
                            std::cout << "\033[31m Cannot find image " <<std::filesystem::path(startImage2).filename().string() << "\033[0m" << std::endl;
                            usage(argv);
                        return 1;
                    }}
					else { // More than two images to open
						std::cout << "\033[31m Too many arguments" << "\033[0m" << std::endl;
						usage(argv);
						return 1;
					}
				} else {
					std::cout << "\033[31m No image-path found" << "\033[0m" << std::endl;
					usage(argv);
					return 1;
				}

			}

			else if (argString == "-scan") { // 3D-scan: Convert with slam6d/fbr/scan_to_panorama
				if (i + 1 < argc) {
					if (startImage.empty() && startScan.empty()) {
						startScan = argv[++i];
                        if(!std::filesystem::exists(startScan)){
                            std::cout << "\033[31m Cannot find scan " <<std::filesystem::path(startScan).filename().string() << "\033[0m" << std::endl;
                            usage(argv);
                        return 1;
                    }
					} else if (startImage2.empty() && startScan2.empty()) {

                        twoImageMode = true;
						startScan2 = argv[++i];
                        if(!std::filesystem::exists(startScan2)){
                            std::cout << "\033[31m Cannot find scan " <<std::filesystem::path(startScan2).filename().string() << "\033[0m" << std::endl;
                            usage(argv);
                        return 1;
                    }
					} else { // More than two images to open
						std::cout << "\033[31m Too many arguments" << "\033[0m" << std::endl;
						usage(argv);
						return 1;
					}
				} else {
					std::cout << "\033[31m No scan-path found" << "\033[0m" << std::endl;
					usage(argv);
					return 1;
				}

			} else if (argString == "-out") { // output directory path for converted scans and coordinates
				if (i + 1 < argc){
					outputDir = argv[++i];
                    if(!std::filesystem::exists(outputDir)){
                        std::cout << "\033[31m Output directory does not exist" << "\033[0m"<< std::endl;
                        usage(argv);
                        return 1;
                    }
                    app.setOutDir(outputDir);}
				else {
					std::cout << "\033[31m No directory path found" << "\033[0m" << std::endl;
					usage(argv);
					return 1;
				}
			} else if (argString == "-help"){
                usage(argv);
                return 1;
            } else{
				std::cout << "\033[31m Invalid argument." << "\033[0m" << std::endl;
				usage(argv);
				return 1;
			}
		}
		// convert scans to 2D-image
		if (!startScan.empty()) {
            
			startImage = app.Create_Panorama(startScan);
            if(startImage.empty()){
                std::cout << "\033[31m" << app.m_convertErrorMessage << "\033[0m" << std::endl;
                usage(argv);
                return 1;
            }}

		
        if (!startScan2.empty()) {
            startImage2 = app.Create_Panorama(startScan2);
            if(startImage2.empty()){
                std::cout << "\033[31m" << app.m_convertErrorMessage << "\033[0m" << std::endl;
                usage(argv);
                return 1;
            }
		}
	}

	// start app
    if (twoImageMode) {
        app.InitTwoImages(startImage, startImage2);
    } else {
        app.Init(startImage);
    }

	while (!glfwWindowShouldClose(window) && !app.ShouldClose()) {
		glfwPollEvents();
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		//
		app.Update();

		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		glfwSwapBuffers(window);
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
