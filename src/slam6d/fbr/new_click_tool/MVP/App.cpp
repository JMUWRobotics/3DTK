#include "App.h"
#include "imgui.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm> // Für std::min

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace fs = std::filesystem;

App::App() {}
App::~App() {
    if (m_texture) glDeleteTextures(1, &m_texture);
}

void App::Init(const std::string& initialImagePath) {
    if (!fs::exists("Koordinaten")) {
        fs::create_directory("Koordinaten");
    }

    if (!initialImagePath.empty()) {
        LoadWorkspace(initialImagePath);
    }
}

void App::LoadWorkspace(const std::string& imagePath) {
    m_currentImagePath = imagePath;
    strncpy(m_imageInputBuffer, imagePath.c_str(), sizeof(m_imageInputBuffer));

    m_imageLoaded = LoadTexture(imagePath);
    
    m_points.clear();

    if (m_imageLoaded) {
        fs::path p(imagePath);
        m_currentTxtPath = "Koordinaten/" + p.stem().string() + "_koordinaten.txt";
        LoadPointsFromFile();
        
        // das Bild im ersten Frame zentrieren und vollständig darstellen 
        m_needsFit = true; 
    }
}

bool App::LoadTexture(const std::string& filename) {
    if (m_texture) glDeleteTextures(1, &m_texture);

    int channels;
    unsigned char* data = stbi_load(filename.c_str(), &m_imgWidth, &m_imgHeight, &channels, 4);
    if (!data) return false;

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
        Point p;
        if (std::getline(ss, val, ',')) p.x = std::stof(val);
        if (std::getline(ss, val, ',')) p.y = std::stof(val);
        m_points.push_back(p);
    }
}

void App::SavePointsToFile() {
    if (m_currentTxtPath.empty()) return;
    std::ofstream file(m_currentTxtPath, std::ios::trunc);
    for (const auto& p : m_points) {
        file << p.x << "," << p.y << "\n";
    }
}

void App::Update() {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();

    // 1. Vollbild-Fenster für das Bild
    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(viewport->WorkSize);
    ImGui::Begin("Workspace", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoBringToFrontOnFocus);
    
    if (m_imageLoaded) {
        ImGuiIO& io = ImGui::GetIO();
        
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

        // Zoom (Mausrad)
        if (ImGui::IsWindowHovered() && io.MouseWheel != 0.0f) {
            // Zoom in die Mitte des Bildschirms justieren
            float zoomFactor = 1.0f + (io.MouseWheel * 0.1f);
            m_zoom *= zoomFactor;
            
            // Angepasste Limits für extrem große Bilder
            if (m_zoom < 0.01f) m_zoom = 0.01f; 
            if (m_zoom > 20.0f) m_zoom = 20.0f;
        }

        ImVec2 p_min = ImVec2(viewport->WorkPos.x + m_panX, viewport->WorkPos.y + m_panY);
        ImVec2 p_max = ImVec2(p_min.x + m_imgWidth * m_zoom, p_min.y + m_imgHeight * m_zoom);

        ImGui::GetWindowDrawList()->AddImage((void*)(intptr_t)m_texture, p_min, p_max);

        if (m_showPoints) {
            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            for (const auto& p : m_points) {
                ImVec2 center = ImVec2(p_min.x + p.x * m_zoom, p_min.y + p.y * m_zoom);
                draw_list->AddCircleFilled(center, 5.0f * m_zoom, IM_COL32(255, 0, 0, 255));
                draw_list->AddCircle(center, 6.0f * m_zoom, IM_COL32(255, 255, 255, 200));
            }
        }

        if (m_selectionMode && ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !ImGui::IsAnyItemHovered()) {
            ImVec2 mousePos = ImGui::GetMousePos();
            float pixelX = (mousePos.x - p_min.x) / m_zoom;
            float pixelY = (mousePos.y - p_min.y) / m_zoom;

            if (pixelX >= 0 && pixelX <= m_imgWidth && pixelY >= 0 && pixelY <= m_imgHeight) {
                m_points.push_back({pixelX, pixelY});
                SavePointsToFile(); 
            }
        }
    } else {
        ImGui::TextColored(ImVec4(1, 0, 0, 1), "Kein Bild geladen oder Bildpfad ungueltig.");
    }
    ImGui::End();

    // Schwebendes UI-Fenster
    ImGui::Begin("Steuerung & Setup", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    
    ImGui::Text("Bild laden:");
    ImGui::InputText("##imagepath", m_imageInputBuffer, sizeof(m_imageInputBuffer));
    ImGui::SameLine();
    if (ImGui::Button("Laden")) {
        LoadWorkspace(m_imageInputBuffer);
    }

    ImGui::Separator();
    
    ImGui::Checkbox("Auswahl-Modus", &m_selectionMode);
    ImGui::SameLine();
    ImGui::Checkbox("Punkte anzeigen", &m_showPoints);
    
    // Button um manuell das Bild wieder zentriert und passend zu machen
    if (ImGui::Button("Ansicht zentrieren & anpassen", ImVec2(-1, 0))) {
        m_needsFit = true;
    }

    ImGui::Separator();
    ImGui::Text("Steuerung:");
    ImGui::Text("- Rechter Mausklick (halten) = Verschieben");
    ImGui::Text("- Mausrad = Zoomen");
    
    ImGui::Separator();
    ImGui::Text("Punkte: %d", (int)m_points.size());
    if (!m_currentTxtPath.empty()) {
        ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.5f, 1.0f), "Ort: %s", m_currentTxtPath.c_str());
    }

    if (ImGui::Button("Rückgängig")) {
        if (!m_points.empty()) {
            m_points.pop_back();
            SavePointsToFile();
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Alle löschen")) {
        m_points.clear();
        SavePointsToFile();
    }

    ImGui::Separator();
    if (ImGui::Button("Speichern & Beenden", ImVec2(-1, 30))) { 
        SavePointsToFile();
        m_shouldClose = true;
    }

    ImGui::End();
}
