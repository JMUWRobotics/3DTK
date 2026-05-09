#pragma once
#include <string>
#include <vector>
#include <glad/glad.h>

struct Point {
    float x, y;
};

class App {
public:
    App();
    ~App();

    void Init(const std::string& initialImagePath);
    void Update();
    bool ShouldClose() const { return m_shouldClose; }

private:
    void LoadWorkspace(const std::string& imagePath);
    void SavePointsToFile();
    void LoadPointsFromFile();
    bool LoadTexture(const std::string& filename);

    // Status-Variablen
    bool m_shouldClose = false;
    bool m_selectionMode = false;
    bool m_showPoints = true;
    bool m_needsFit = false; // Sagt der App, dass sie im nächsten Frame zoomen muss
    
    // Bild-Daten
    std::string m_currentImagePath;
    std::string m_currentTxtPath;
    GLuint m_texture = 0;
    int m_imgWidth = 0;
    int m_imgHeight = 0;
    bool m_imageLoaded = false;

    // View-Daten
    float m_zoom = 1.0f;
    float m_panX = 0.0f;
    float m_panY = 0.0f;

    // Gespeicherte Punkte
    std::vector<Point> m_points;

    // Puffer für die UI-Texteingabe
    char m_imageInputBuffer[256] = "";
};
