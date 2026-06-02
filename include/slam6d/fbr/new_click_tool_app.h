#pragma once
#include <string>
#include <vector>
// #include <glad/glad.h>

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif

struct ClickPoint {
	float x, y;
};

struct Correspondence {
    ClickPoint first;
    ClickPoint second;
};

class App
{
      public:
	App();
	~App();

    std::string Create_Panorama(const std::string& startScan);
    void setOutDir(std::string outputDir);
	void Init(const std::string &initialImagePath);
    void InitTwoImages(const std::string& firstImagePath, const std::string& secondImagePath);
	void Update();
	bool ShouldClose() const { return m_shouldClose; }

      private:
	void LoadWorkspace(const std::string &imagePath);
    void LoadTwoImageWorkspace(const std::string& firstImagePath,
                           const std::string& secondImagePath);
	void ResetWorkspace();


	void SavePointsToFile();
	void LoadPointsFromFile();
	bool LoadTexture(const std::string &filename);
    bool LoadTextureFromMemory(const unsigned char* data, int width, int height); //two image


	// Status-Variablen
    bool m_twoImageMode = false;
	bool m_firstImageIsScan = false;
	bool m_secondImageIsScan = false;
	bool m_shouldClose = false;
	bool m_selectionMode = false;
	bool m_showPoints = true;
	bool m_needsFit = false; // Sagt der App, dass sie im nächsten Frame zoomen muss
	float m_closeButtonSize = 0;
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
	std::vector<ClickPoint> m_points;

        //variables for two image mode:
    int m_firstOffsetX = 0;
    int m_firstOffsetY = 0;
    int m_firstWidth = 0;
    int m_firstHeight = 0;
    int m_secondOffsetX = 0;
    int m_secondOffsetY = 0;
    int m_secondWidth = 0;
    int m_secondHeight = 0;
    std::vector<Correspondence> m_correspondences;
    ClickPoint m_pendingFirstPoint;
    bool m_waitingForSecondPoint = false;

	// Puffer für die UI-Texteingabe
	char m_imageInputBuffer[256] = "";
	char m_imageInputBuffer2[256] = "";


	// OutputDir
	std::string m_outputDir = "";
    char m_outDirBuffer[256] = "";

	public:
		//Error Messages
    std::string m_errorMessage;
	std::string m_inputErrorMessage;
	std::string m_inputErrorMessage2;
	std::string m_convertErrorMessage;
	std::string m_outputDirErrorMessage = "";

};
