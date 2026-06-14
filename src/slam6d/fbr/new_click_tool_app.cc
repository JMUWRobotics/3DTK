#include "slam6d/fbr/new_click_tool_app.h"
#include "imgui.h"
#include "slam6d/fbr/panorama.h"
#include <algorithm> // Für std::min
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <unistd.h>
#include <sys/wait.h>
#include <GLFW/glfw3.h>
#include <thread>
#include <chrono>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

int run(int argc, char **argv);

namespace fs = std::filesystem;

App::App() {}
App::~App()
{
	if (m_texture)
		glDeleteTextures(1, &m_texture);
}

std::string App::Create_Panorama(const std::string &startScan, const std::string &scanformat)
{
    m_convertErrorMessage = "";
    try{
	std::filesystem::path p(startScan);
	std::string scanDir = p.parent_path().string();
	std::string scanName = p.stem().string();
	std::string scanOutDir = m_outputDir.empty() ? p.parent_path().string() : m_outputDir;
	int scanNr = std::stoi(p.stem().string().substr(4));

    //  command line for scan_to_panorama with normalized range
	std::string command = "scan_to_panorama " + scanDir + " -s " + std::to_string(scanNr) + " -e " +
			      std::to_string(scanNr) + " -f " + scanformat + " -A -a -F PNG -O " + scanOutDir + "\0";
std::cout << command << std::endl;
    // convert command line to string array
	std::stringstream commandStream(command);
	std::vector<std::string> args;
	std::string argTemp;

	while (commandStream >> argTemp) {
		args.push_back(argTemp);
	}

	// vector for scan_to_panorama-input
	std::vector<char *> args_char;
	for (std::string &a : args) {
		args_char.push_back(const_cast<char*>(a.c_str()));
        
	}
    args_char.push_back(nullptr);


	// Use scan_to_panorama in new process to avoid effects of global variables

    pid_t pid = fork();

    if (pid == 0)
    {
	    run(args.size(), args_char.data());
        std::cout << "ende scan to panorama "  << std::endl;

        _exit(0);
    }
    // avoids window-freeze messages while converting a large scan
    int status;
    while(waitpid(pid, &status, WNOHANG) == 0){
        glfwPollEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

	// check if conversion worked
	std::string genImName = scanOutDir + "/" + scanName + "_EQUIRECTANGULAR_3600x1000_NormalizedRange.png";
	if (!std::filesystem::exists(genImName)) {

        m_convertErrorMessage = "Failed generating panorama from " + startScan;
		return "";
	} else {

		std::cout << "Panorama created in " << scanOutDir << std::endl << std::endl;
		return genImName;
	}} catch(const std::exception& e){

        m_convertErrorMessage = "Failed generating panorama from " + startScan;
        return "";
    }
}

void App::setOutDir(std::string outputDir) { m_outputDir = outputDir; }

void App::Init(const std::string &initialImagePath)
{
   
	if (!initialImagePath.empty()) {
		LoadWorkspace(initialImagePath);
	}
}

void App::LoadWorkspace(const std::string &imagePath)
{
	m_currentImagePath = imagePath;
	strncpy(m_imageInputBuffer, imagePath.c_str(), sizeof(m_imageInputBuffer)-1);

	m_imageLoaded = LoadTexture(imagePath);
	m_points.clear();

	if (m_imageLoaded) {
		fs::path p(imagePath);

        m_outputDir = m_outputDir.empty() ? p.parent_path().string() : m_outputDir;


		m_currentTxtPath = m_outputDir + "/Koordinaten/" + p.stem().string() + "_koordinaten.txt";
		LoadPointsFromFile();
		// das Bild im ersten Frame zentrieren und vollständig darstellen
		m_needsFit = true;
	}
}

void App::ResetWorkspace(){
    m_imageLoaded = false;
    m_twoImageMode = false;
    m_firstImageIsScan = false;
    m_secondImageIsScan = false;
    m_zoom = 1.0f;
    m_panX = 0.0f;
    m_panY = 0.0f;
    m_points.clear();
    m_correspondences.clear();
    m_errorMessage.clear();
    m_inputErrorMessage.clear();
    m_inputErrorMessage2.clear();
    m_waitingForSecondPoint = false;

    m_imageInputBuffer[0] = '\0';
    m_imageInputBuffer2[0] = '\0';
    m_currentImagePath.clear();
    m_currentTxtPath.clear();

    if(m_texture != 0){
        glDeleteTextures(1, &m_texture);
        m_texture = 0;}

}

bool App::LoadTexture(const std::string &filename)
{
	if (m_texture)
		glDeleteTextures(1, &m_texture);

	int channels;
	unsigned char *data = stbi_load(filename.c_str(), &m_imgWidth, &m_imgHeight, &channels, 4);
	if (!data)
		return false;

	glGenTextures(1, &m_texture);
	glBindTexture(GL_TEXTURE_2D, m_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_imgWidth, m_imgHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
	stbi_image_free(data);
	return true;
}

void App::LoadPointsFromFile() {
    if (!fs::exists(m_currentTxtPath)) return; 

    std::ifstream file(m_currentTxtPath);
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string val;
        if (m_twoImageMode) {
            Correspondence c;

            if (!std::getline(ss, val, ',')) continue;
            c.first.x = std::stof(val);

            if (!std::getline(ss, val, ',')) continue;
            c.first.y = std::stof(val);

            if (!std::getline(ss, val, ',')) continue;
            c.second.x = std::stof(val);

            if (!std::getline(ss, val, ',')) continue;
            c.second.y = std::stof(val);

            m_correspondences.push_back(c);
        } else {
            ClickPoint p;

            if (!std::getline(ss, val, ',')) continue;
            p.x = std::stof(val);

            if (!std::getline(ss, val, ',')) continue;
            p.y = std::stof(val);

            m_points.push_back(p);
        }
    }
}

void App::SavePointsToFile() {
    if (m_currentTxtPath.empty()) return;
    //deletes all previous information and saves the new
    std::ofstream file(m_currentTxtPath, std::ios::trunc);

    if (m_twoImageMode) {
        for (const auto& c : m_correspondences) {
            file << c.first.x << "," << c.first.y << ","
                 << c.second.x << "," << c.second.y << "\n";
        }
    } else {
        for (const auto& p : m_points) {
            file << p.x << "," << p.y << "\n";
        }
    }
}

void App::Update()
{
	const ImGuiViewport *viewport = ImGui::GetMainViewport();

	// 1. Vollbild-Fenster für das Bild
	ImGui::SetNextWindowPos(viewport->WorkPos);
	ImGui::SetNextWindowSize(viewport->WorkSize);
	ImGui::Begin("Workspace", nullptr,
		     ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoBackground |
			 ImGuiWindowFlags_NoBringToFrontOnFocus);
	if (m_imageLoaded) {
		ImGuiIO &io = ImGui::GetIO();
		// --- Auto-Fit Logik ---
		if (m_needsFit) {
			// Berechne, wie stark wir auf X und Y zoomen müssten, damit es passt
			float zoomX = viewport->WorkSize.x / (float)m_imgWidth;
			float zoomY = viewport->WorkSize.y / (float)m_imgHeight;
			// den kleineren Zoom nehmen, damit nichts abgeschnitten wird
			m_zoom = std::min(zoomX, zoomY);
			// Mache es ein bisschen kleiner, damit es nicht exakt am Rand klebt
			m_zoom *= 0.95f;

			// Zentriere das Bild
			m_panX = (viewport->WorkSize.x - (m_imgWidth * m_zoom)) * 0.5f;
			m_panY = (viewport->WorkSize.y - (m_imgHeight * m_zoom)) * 0.5f;
			m_needsFit = false;
		}
		// ---------------------------

		// Verschieben
		if (ImGui::IsWindowHovered() && ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
			m_panX += io.MouseDelta.x;
			m_panY += io.MouseDelta.y;
		}

		// Zoom to mouse position
		if (ImGui::IsWindowHovered() && io.MouseWheel != 0.0f) {
			// current mouse / pixel-position
			ImVec2 mousePosition = ImGui::GetMousePos();
			float curX = (mousePosition.x - (viewport->WorkPos.x + m_panX)) / m_zoom;
			float curY = (mousePosition.y - (viewport->WorkPos.y + m_panY)) / m_zoom;
			float zoomFactor = 1.0f + (io.MouseWheel * 0.1f);
			m_zoom *= zoomFactor;
			// Angepasste Limits für extrem große Bilder
			if (m_zoom < 0.01f)
				m_zoom = 0.01f;
			if (m_zoom > 20.0f)
				m_zoom = 20.0f;
			// fix pan-position:
			m_panX = mousePosition.x - viewport->WorkPos.x - curX * m_zoom;
			m_panY = mousePosition.y - viewport->WorkPos.y - curY * m_zoom;
		}

		ImVec2 p_min = ImVec2(viewport->WorkPos.x + m_panX, viewport->WorkPos.y + m_panY);
		ImVec2 p_max = ImVec2(p_min.x + m_imgWidth * m_zoom, p_min.y + m_imgHeight * m_zoom);

        ImGui::GetWindowDrawList()->AddImage((void*)(intptr_t)m_texture, p_min, p_max);

        if (m_showPoints) {
             ImDrawList* draw_list = ImGui::GetWindowDrawList();

            if (m_twoImageMode) {
                for (const auto& c : m_correspondences) {
                    ImVec2 firstCenter = ImVec2(
                        p_min.x + (m_firstOffsetX + c.first.x) * m_zoom,
                        p_min.y + (m_firstOffsetY + c.first.y) * m_zoom
                    );

                    ImVec2 secondCenter = ImVec2(
                        p_min.x + (m_secondOffsetX + c.second.x) * m_zoom,
                        p_min.y + (m_secondOffsetY + c.second.y) * m_zoom
                    );

                    draw_list->AddLine(firstCenter, secondCenter, IM_COL32(255, 220, 0, 255), 2.0f);
                    draw_list->AddCircleFilled(firstCenter, 5.0f * m_zoom, IM_COL32(255, 0, 0, 255));
                    draw_list->AddCircle(firstCenter, 6.0f * m_zoom, IM_COL32(255, 255, 255, 200));

                    draw_list->AddCircleFilled(secondCenter, 5.0f * m_zoom, IM_COL32(0, 180, 255, 255));
                    draw_list->AddCircle(secondCenter, 6.0f * m_zoom, IM_COL32(255, 255, 255, 200));
                }
            } else {
                for (const auto& p : m_points) {
                    ImVec2 center = ImVec2(p_min.x + p.x * m_zoom, p_min.y + p.y * m_zoom);
                    draw_list->AddCircleFilled(center, 5.0f * m_zoom, IM_COL32(255, 0, 0, 255));
                    draw_list->AddCircle(center, 6.0f * m_zoom, IM_COL32(255, 255, 255, 200));
                }
            }
        }

        if (m_selectionMode && ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !ImGui::IsAnyItemHovered()) {
            ImVec2 mousePos = ImGui::GetMousePos();
            float pixelX = (mousePos.x - p_min.x) / m_zoom;
            float pixelY = (mousePos.y - p_min.y) / m_zoom;

            if (m_twoImageMode) {
                bool insideFirst =
                    pixelX >= m_firstOffsetX &&
                    pixelX <= m_firstOffsetX + m_firstWidth &&
                    pixelY >= m_firstOffsetY &&
                    pixelY <= m_firstOffsetY + m_firstHeight;

                bool insideSecond =
                    pixelX >= m_secondOffsetX &&
                    pixelX <= m_secondOffsetX + m_secondWidth &&
                    pixelY >= m_secondOffsetY &&
                    pixelY <= m_secondOffsetY + m_secondHeight;

                if (!m_waitingForSecondPoint) {
                    if (insideFirst) {
                        m_pendingFirstPoint = {
                            pixelX - m_firstOffsetX,
                            pixelY - m_firstOffsetY
                        };
                        m_waitingForSecondPoint = true;
                        m_errorMessage.clear();
                    } else {
                        m_errorMessage = "Select point in first image.";
                    }
                } else {
                    if (insideSecond) {
                        ClickPoint secondPoint = {
                            pixelX - m_secondOffsetX,
                            pixelY - m_secondOffsetY
                        };

                        m_correspondences.push_back({
                            m_pendingFirstPoint,
                            secondPoint
                        });
                        SavePointsToFile();
                        m_waitingForSecondPoint = false;
                        m_errorMessage.clear();
                    } else {
                        m_errorMessage = "Select a point in second image.";
                    }
                }
            } else {
                if (pixelX >= 0 && pixelX <= m_imgWidth && pixelY >= 0 && pixelY <= m_imgHeight) {
                    m_points.push_back({pixelX, pixelY});
                    SavePointsToFile();
                }
            }
        }
    } else {
        ImGui::TextColored(ImVec4(1, 0, 0, 1), "No image loaded");
    }
    ImGui::End();

	// Schwebendes UI-Fenster
	ImGui::Begin("Control & Setup", nullptr, ImGuiWindowFlags_AlwaysAutoResize);


//set output dir
        //Disable set output and load-section when images successfully created
    ImGui::BeginDisabled(m_imageLoaded);

if(ImGui::Button("Set out-dir")){
        ImGui::OpenPopup("Set output directory");}
        if(ImGui::BeginPopupModal("Set output directory", nullptr, ImGuiWindowFlags_AlwaysAutoResize)){
        ImGui::Text("output directory:");
        ImGui::InputText("##SaveDirectory", m_outDirBuffer, sizeof(m_outDirBuffer));
        if(ImGui::Button("Save")){
            if(fs::exists(m_outDirBuffer)){
                if(fs::is_directory(m_outDirBuffer)){
            m_outputDirErrorMessage = "";
            m_outputDir = m_outDirBuffer;
            ImGui::CloseCurrentPopup();}
            else{ if(fs::exists(m_outDirBuffer))
                m_outputDirErrorMessage = "Not a directory";
            }
        }
            else{
                m_outputDirErrorMessage = "Path does not exist";
            }}
        ImGui::SameLine();
        if(ImGui::Button("Cancel")){
        ImGui::CloseCurrentPopup();  
        }
        ImGui::TextColored(ImVec4(1.0f, 0.2f, 0.2f, 1.0f), "%s", m_outputDirErrorMessage.c_str());

        ImGui::EndPopup();
    }
        //help-function
    ImGui::SameLine();
    ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - 20);
   if(ImGui::Button("?")){
        ImGui::OpenPopup("Help");
    }
    if(ImGui::BeginPopupModal("Help", nullptr, ImGuiWindowFlags_AlwaysAutoResize)){
	    ImGui::Text("Control:");
	    ImGui::Text("- Right mouseclick (hold) = slide picture");
	    ImGui::Text("- Mousewheel = zoom");
        ImGui::Separator();

/////////////////////////////// TODO: Anweisungen/Hilfestellungen einfügen



        ImGui::Separator();

    ImGui::SetCursorPosX((ImGui::GetWindowContentRegionMax().x - m_closeButtonSize)/2);
    if(ImGui::Button("Close")){
        ImGui::CloseCurrentPopup();
    }
        m_closeButtonSize = ImGui::GetItemRectSize().x;

        ImGui::EndPopup();
    }

       //Load-section:


	ImGui::InputTextWithHint("##imagepath", m_twoImageMode ? "path to first image" : "path to image", m_imageInputBuffer, sizeof(m_imageInputBuffer));
    ImGui::SameLine();
    ImGui::Checkbox("Scan", &m_firstImageIsScan);
    if(m_firstImageIsScan){
        ImGui::SameLine();
    if (ImGui::BeginCombo("Format", m_current_item)){
        for(int i = 0; i < IM_ARRAYSIZE(m_formatitems); i++){
            bool is_selected = m_current_item == m_formatitems[i];
            if(ImGui::Selectable(m_formatitems[i], is_selected)) m_current_item = m_formatitems[i];
            if(is_selected) ImGui::SetItemDefaultFocus();

        }
        ImGui::EndCombo();

    }}
    
    if(!m_inputErrorMessage.empty())ImGui::TextColored(ImVec4(1, 0, 0, 1), "%s", m_inputErrorMessage.c_str());

        //one or two images mode
    ImGui::Checkbox("two-image mode", &m_twoImageMode);
    if(m_twoImageMode){
        ImGui::InputTextWithHint("##imagepath2", "path to second image", m_imageInputBuffer2, sizeof(m_imageInputBuffer2));
        ImGui::SameLine();
    	ImGui::Checkbox("Scan##2", &m_secondImageIsScan);
    if(m_secondImageIsScan){
        ImGui::SameLine();
    if (ImGui::BeginCombo("Format##2", m_current_item2)){
        for(int i = 0; i < IM_ARRAYSIZE(m_formatitems); i++){
            bool is_selected = m_current_item2 == m_formatitems[i];
            if(ImGui::Selectable(m_formatitems[i], is_selected)) m_current_item2 = m_formatitems[i];
            if(is_selected) ImGui::SetItemDefaultFocus();

        }
        ImGui::EndCombo();

    }}
    if(!m_inputErrorMessage2.empty())ImGui::TextColored(ImVec4(1, 0, 0, 1), "%s", m_inputErrorMessage2.c_str());

    }
    bool error = false;
    //disable load-button if too few arguments given 
        ImGui::BeginDisabled(strlen(m_imageInputBuffer) == 0 ||( m_twoImageMode && strlen(m_imageInputBuffer2) == 0) || m_conversion_ongoing);
        if (ImGui::Button(m_twoImageMode ? "Load images" : "Load image", ImVec2(-1, 30))) {
        // if scan: convert
        m_inputErrorMessage = "";
        
        if(m_firstImageIsScan){
            std::string imInputBufferString;
            if(!fs::exists(m_imageInputBuffer)){
                error = true;
                m_inputErrorMessage = "Cannot find scan " + std::filesystem::path(m_imageInputBuffer).filename().string();}
            else{ imInputBufferString = Create_Panorama(m_imageInputBuffer, m_current_item);
            if(imInputBufferString == ""){
                m_inputErrorMessage =  m_convertErrorMessage;
                error = true;
	        }
            strncpy(m_imageInputBuffer, imInputBufferString.c_str(), sizeof(m_imageInputBuffer)-1);
        }
        } else {    //first image is image: test if exists
                if(!fs::exists(m_imageInputBuffer)){
                error = true;
                m_inputErrorMessage = "Cannot find image " + std::filesystem::path(m_imageInputBuffer).filename().string();}
        }
        if(m_twoImageMode){
            m_inputErrorMessage2 = "";
            if(m_secondImageIsScan){
            std::string imInputBuffer2String;
            if(!fs::exists(m_imageInputBuffer2)){
                error = true;
                m_inputErrorMessage2 = "Cannot find scan " + std::filesystem::path(m_imageInputBuffer2).filename().string();}
            else{imInputBuffer2String = Create_Panorama(m_imageInputBuffer2, m_current_item2);
            if(imInputBuffer2String == ""){
                m_inputErrorMessage2 =  m_convertErrorMessage;
                error = true;}
	        
            strncpy(m_imageInputBuffer2, imInputBuffer2String.c_str(), sizeof(m_imageInputBuffer2)-1);
            } }else {    //second image is image: test if exists
                if(!fs::exists(m_imageInputBuffer2)){
                error = true;
                m_inputErrorMessage2 = "Cannot find image " + std::filesystem::path(m_imageInputBuffer2).filename().string();}
        }}

        
        if(!error){
        if(m_twoImageMode) LoadTwoImageWorkspace(m_imageInputBuffer, m_imageInputBuffer2);
        else LoadWorkspace(m_imageInputBuffer);
    	}
        if(!m_imageLoaded) m_errorMessage = "Failed loading workspace"; 
}
    ImGui::EndDisabled();

    ImGui::EndDisabled();

    //Reset Workspace
    ImGui::BeginDisabled(!m_imageLoaded);

    if(ImGui::Button("Reset workspace", ImVec2(-1, 0))){
        ImGui::OpenPopup("Save points?");}
        if(ImGui::BeginPopupModal("Save points?", nullptr, ImGuiWindowFlags_AlwaysAutoResize)){
        ImGui::Text("Save points before reset?");
        if(ImGui::Button("Save")){
        SavePointsToFile();    
        ResetWorkspace();
        ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if(ImGui::Button("Don't save")){
        ResetWorkspace();
        ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if(ImGui::Button("Cancel")){
        ImGui::CloseCurrentPopup();  
        }
        ImGui::EndPopup();
        }
    

    ImGui::EndDisabled();


	ImGui::Separator();
	ImGui::Checkbox("Selection-mode", &m_selectionMode);
	ImGui::SameLine();
	ImGui::Checkbox("Show points", &m_showPoints);
	// Button um manuell das Bild wieder zentriert und passend zu machen
	if (ImGui::Button("Center view", ImVec2(-1, 0))) {
		m_needsFit = true;
	}


	    ImGui::Separator();
    if (m_twoImageMode) {
        ImGui::Text("Correspondences: %d", (int)m_correspondences.size());

        if (m_selectionMode) {
            if (m_waitingForSecondPoint) {
                ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f),
                                "Choose point in second image.");
            } else {
                ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f),
                                "Choose point in first image.");
            }
        }
    } else {
        ImGui::Text("points: %d", (int)m_points.size());
    }    
    if (!m_currentTxtPath.empty()) {
        ImGui::PushTextWrapPos(400.0f);
        ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.5f, 1.0f), "Location: %s", m_currentTxtPath.c_str());
        ImGui::PopTextWrapPos();
    }

    if (!m_errorMessage.empty()) {
    ImGui::TextColored(ImVec4(1.0f, 0.2f, 0.2f, 1.0f), "%s", m_errorMessage.c_str());
    }
    if(m_imageLoaded){
    if(ImGui::Button("Copy Path")){
        ImGui::SetClipboardText(m_currentTxtPath.c_str());
    }}
    
    if (ImGui::Button("Undo")) {
      if (m_twoImageMode) {
            if (m_waitingForSecondPoint) {
                m_waitingForSecondPoint = false;
                m_errorMessage.clear();
            } else if (!m_correspondences.empty()) {
                m_correspondences.pop_back();
                SavePointsToFile();
            }
        } else {
            if (!m_points.empty()) {
                m_points.pop_back();
                SavePointsToFile();
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("delete all")) {
         if (m_twoImageMode) {
            m_correspondences.clear();
            m_waitingForSecondPoint = false;
            m_errorMessage.clear();
            SavePointsToFile();
        } else {
            m_points.clear();
            SavePointsToFile();
        }
    }

	ImGui::Separator();
	if (ImGui::Button("save & exit", ImVec2(-1, 30))) {
		SavePointsToFile();
		m_shouldClose = true;
	}
	ImGui::End();
}
