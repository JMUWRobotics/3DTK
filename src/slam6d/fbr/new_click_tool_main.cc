// #define GL_SILENCE_DEPRECATION

// #include <glad/glad.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <filesystem>
#include <iostream>

#include "slam6d/fbr/new_click_tool_app.h"

//usage function to print in the console window
void usage(char **argv)
{
	std::cout << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << std::endl;
	std::cout << "-help \t\t\t\t\t\t\t\t\t Option list" << std::endl;
    std::cout << std::endl;
	std::cout << "<No Arguments> \t\t\t\t\t\t\t\t Open empty window" << std::endl;
	std::cout << "<input format> <path> \t\t\t\t\t\t\t Open one image, default settings" << std::endl;
	std::cout << "<input format> <path> <input format> <path> \t\t\t\t Open two images, default settings" << std::endl;
	std::cout << "<input format> <path> <out dir> \t\t\t\t\t Open one image, default settings, manually selected output directory" << std::endl;

	std::cout << std::endl;
	std::cout << "-im <image path> \t\t\t\t\t\t\t No further settings if loading 2D-image" << std::endl; 
	std::cout << "-scan <scan path> -f <scan format>\t\t\t\t\t Open one image based on 3D-scan, manually selected format, default directory and conversion mode" << std::endl; 
	std::cout << "-im <image path> -scan <scan path> -f <scan format> <conversion mode> \t Open image and 3D-scan image, manually selected format and conversion mode, default directory " << std::endl; 
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "input format:" << std::endl;
	std::cout << "-im \t\t\t\t\t\t 2D-image Format" << std::endl;
	std::cout << "-scan \t\t\t\t\t\t 3D-scan Format, internally converted with scan_to_panorama" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "conversion mode:" << std::endl;
	std::cout << "-A \t\t\t\t\t\t Range" << std::endl;
	std::cout << "-R \t\t\t\t\t\t Reflectance" << std::endl;
	std::cout << "-a \t\t\t\t\t\t Normalized Range" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "Possible scan formats:" << "\t\t\t\t uos, uosc, uos_map, uos_rgb, uos_frames, uos_map_frames, old, rts, rts_map, ifp, riegl_txt,\n \t\t\t\t\t\t riegl_rgb, riegl_bin, zahn, ply, wrl, xyz, xyzc, zuf, iais, front, x3d, rxp, ais" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "Default scan format\t\t\t\t uos" << std::endl;
	std::cout << "Default conversion mode\t\t\t\t Normalized Range" << std::endl;
	std::cout << "Default output directory\t\t\t directory of first loaded image/scan" << std::endl;
    std::cout << "Converted scans\t\t\t\t\t output directory or directory of scan" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Example:\t\t\t\t\t -im /pictures/myPicture.png -scan /scans/scan001.3D -f -R uosr -out /users/desktop/click_tool \n"
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
	
	GLFWmonitor* primary = glfwGetPrimaryMonitor();
	//const GLFWvidmode* mode = glfwGetVideoMode(primary);
	int xpos, ypos, width, height;
glfwGetMonitorWorkarea(primary, &xpos, &ypos, &width, &height);


	GLFWwindow *window = glfwCreateWindow(width,height, "New click tool", NULL, NULL);
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
	std::string scan_format = app.m_formatitems[0]; //Default scanformat = uos
	std::string scan_format2 = app.m_formatitems[0];
    bool twoImageMode = false;
	bool loadFromTerminal = false;

	if (argc > 1) { // if no picture is loaded start app without initial images
		loadFromTerminal = true;
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

			}else if (argString == "-f"){ //Choose scan format if different from ous-Format
				if (i + 1 < argc) {
					std::string format = argv[++i];

					// check if format is valid
				bool format_valid = false;

				for(int i = 0; i < IM_ARRAYSIZE(app.m_formatitems); i++){
					if(format == app.m_formatitems[i]) {format_valid = true;
					break;}
				}
					if(!format_valid){
					std::cout << "\033[31m Format " + format + " is not a valid format" << "\033[0m" << std::endl;
					usage(argv);
					return 1;
				}
				if (startImage.empty() && startScan2.empty()) scan_format = format;
				else if (startImage2.empty() && !startScan2.empty()) scan_format2 = format;
				else{
					std::cout << "\033[31m Too many arguments" << "\033[0m" << std::endl;
						usage(argv);
						return 1;
				}
				


			}else {
					std::cout << "\033[31m Insert valid format" << "\033[0m" << std::endl;
					usage(argv);
					return 1;
				}

			} else if(argString == "-R" ||argString == "-A"|| argString == "-a"){
				if (startImage.empty() && startScan2.empty()) app.m_Conversion = argString;
				else if (startImage2.empty() && !startScan2.empty()) app.m_Conversion2 = argString;
				else{
					std::cout << "\033[31m Too many arguments" << "\033[0m" << std::endl;
						usage(argv);
						return 1;
				}
			}
			
			
			
			
			else if (argString == "-out") { // output directory path for converted scans and coordinates
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
            
			startImage = app.Create_Panorama(startScan, scan_format, app.m_Conversion);
            if(startImage.empty()){
                std::cout << "\033[31m" << app.m_convertErrorMessage << "\033[0m" << std::endl;
                usage(argv);
                return 1;
            }}

		
        if (!startScan2.empty()) {
            startImage2 = app.Create_Panorama(startScan2, scan_format2, app.m_Conversion2);
            if(startImage2.empty()){
                std::cout << "\033[31m" << app.m_convertErrorMessage << "\033[0m" << std::endl;
                usage(argv);
                return 1;
            }
		}
	}



	while (!glfwWindowShouldClose(window) && !app.ShouldClose()) {
		glfwPollEvents();
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		//
		app.Update();

		if(loadFromTerminal){
		// start app
    	if (twoImageMode) {
        	app.InitTwoImages(startImage, startImage2);
    	} else {
        	app.Init(startImage);
    	}

		loadFromTerminal = false;
		}


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
