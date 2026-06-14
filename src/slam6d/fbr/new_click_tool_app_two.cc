#include "slam6d/fbr/new_click_tool_app.h"
#include "imgui.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm> // Für std::min

#include "stb_image.h"

namespace fs = std::filesystem;

/*
The purpose of this file is having some functions only used
for the two-image mode, but there are still some things of that mode
that can be found on App.cpp in functions: LoadPointsFromFile,
SavePointsToFile and Update
*/

void App::InitTwoImages(const std::string& firstImagePath, const std::string& secondImagePath) {
    LoadTwoImageWorkspace(firstImagePath, secondImagePath);
}

void App::LoadTwoImageWorkspace(const std::string& firstImagePath, //Luis
                                const std::string& secondImagePath) {
    m_twoImageMode = true;
    m_currentImagePath = firstImagePath + " | " + secondImagePath;
    strncpy(m_imageInputBuffer, m_currentImagePath.c_str(), sizeof(m_imageInputBuffer)-1);
    m_errorMessage.clear();
    m_waitingForSecondPoint = false;                            
    int firstWidth = 0;
    int firstHeight = 0;
    int firstChannels = 0;

    int secondWidth = 0;
    int secondHeight = 0;
    int secondChannels = 0;

    //load images
    //channels tell us the real number of channels, if RGBA, 4, eventhough we always put 4
    unsigned char* firstData = stbi_load(firstImagePath.c_str(), &firstWidth, &firstHeight, &firstChannels, 4);
    unsigned char* secondData = stbi_load(secondImagePath.c_str(), &secondWidth, &secondHeight, &secondChannels, 4);

    if (!firstData || !secondData) {
        if (firstData) stbi_image_free(firstData);
        if (secondData) stbi_image_free(secondData);

        m_imageLoaded = false;
        return;
    }

        //Create Path to Coordinates.Directory if not yet determined
        fs::path p(firstImagePath);

        m_outputDir = m_outputDir.empty() ? p.parent_path().string() : m_outputDir;

		m_currentTxtPath = m_outputDir + "/Koordinaten/" + p.stem().string() + "_koordinaten.txt";


    //this establishes where each image is positionated
    const int gap = 20;

        //screen size
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    float width = viewport -> WorkSize.x;
    float height = viewport -> WorkSize.y;
    float screenRelation = width/height;
        //imafe size
    float horImageRelation = (float)(firstWidth + gap + secondWidth) / std::max(firstHeight, secondHeight);
    float vertImageRelation = (float)std::max(firstWidth, secondWidth) / (firstHeight + gap +  secondHeight);
        //alignment horizontal/vertical
    bool horizontal = std::abs(screenRelation - horImageRelation) < std::abs(screenRelation - vertImageRelation);

    int combinedWidth = 0;
    int combinedHeight = 0;

    int firstOffsetX = 0;
    int firstOffsetY = 0;
    int secondOffsetX = 0;
    int secondOffsetY = 0;

    if (horizontal) {
        combinedWidth = firstWidth + gap + secondWidth;
        combinedHeight = std::max(firstHeight, secondHeight);

        firstOffsetX = 0;
        firstOffsetY = (combinedHeight - firstHeight) / 2;

        secondOffsetX = firstWidth + gap;
        secondOffsetY = (combinedHeight - secondHeight) / 2;
    } else {
        combinedWidth = std::max(firstWidth, secondWidth);
        combinedHeight = firstHeight + gap + secondHeight;

        firstOffsetX = (combinedWidth - firstWidth) / 2;
        firstOffsetY = 0;

        secondOffsetX = (combinedWidth - secondWidth) / 2;
        secondOffsetY = firstHeight + gap;
    }
    //asignation into de class
    m_firstOffsetX = firstOffsetX;
    m_firstOffsetY = firstOffsetY;
    m_firstWidth = firstWidth;
    m_firstHeight = firstHeight;

    m_secondOffsetX = secondOffsetX;
    m_secondOffsetY = secondOffsetY;
    m_secondWidth = secondWidth;
    m_secondHeight = secondHeight;

    std::vector<unsigned char> combined(combinedWidth * combinedHeight * 4, 0);

    // gray background
    for (int y = 0; y < combinedHeight; ++y) {
        for (int x = 0; x < combinedWidth; ++x) {
            int idx = (y * combinedWidth + x) * 4;
            combined[idx + 0] = 35;
            combined[idx + 1] = 35;
            combined[idx + 2] = 35;
            combined[idx + 3] = 255;
        }
    }

    //here we fullfil combined
    auto copyImage = [](const unsigned char* src,
                        int srcWidth,
                        int srcHeight,
                        std::vector<unsigned char>& dst,
                        int dstWidth,
                        int offsetX,
                        int offsetY) {
        for (int y = 0; y < srcHeight; ++y) {
            for (int x = 0; x < srcWidth; ++x) {
                int srcIdx = (y * srcWidth + x) * 4;
                int dstIdx = ((y + offsetY) * dstWidth + (x + offsetX)) * 4;

                dst[dstIdx + 0] = src[srcIdx + 0];
                dst[dstIdx + 1] = src[srcIdx + 1];
                dst[dstIdx + 2] = src[srcIdx + 2];
                dst[dstIdx + 3] = src[srcIdx + 3];
            }
        }
    };

    copyImage(firstData, firstWidth, firstHeight, combined, combinedWidth, firstOffsetX, firstOffsetY);
    copyImage(secondData, secondWidth, secondHeight, combined, combinedWidth, secondOffsetX, secondOffsetY);

    m_imageLoaded = LoadTextureFromMemory(combined.data(), combinedWidth, combinedHeight);

    stbi_image_free(firstData);
    stbi_image_free(secondData);

    m_correspondences.clear();
    m_waitingForSecondPoint = false;

    if (m_imageLoaded) {
        fs::path firstPath(firstImagePath);
        fs::path secondPath(secondImagePath);

        m_currentTxtPath = m_outputDir + "/Koordinaten/" +
            firstPath.stem().string() + "_" +
            secondPath.stem().string() +
            "_correspondencies.txt";

        LoadPointsFromFile();
        m_needsFit = true;
    }
}

bool App::LoadTextureFromMemory(const unsigned char* data, int width, int height) { //Luis
    if (m_texture) glDeleteTextures(1, &m_texture);

    if (!data || width <= 0 || height <= 0) {
        return false;
    }

    m_imgWidth = width;
    m_imgHeight = height;

    glGenTextures(1, &m_texture);
    glBindTexture(GL_TEXTURE_2D, m_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_imgWidth, m_imgHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

    return true;
}
