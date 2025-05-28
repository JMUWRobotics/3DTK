#pragma once

#include "calibration/Detector.h"

#include "paruco.hpp"

namespace calibration {

class PArUcoDetector : public Detector {
public:
    PArUcoDetector(float tagExpandScale = 0.5f, const std::string& dictionaryName = "DICT_6X6_250");
    ~PArUcoDetector() = default;
    bool detect(const cv::Mat& image) noexcept(false);
    void writeDetectionsToFile(const std::string& path);
    void readDetectionsFromFile(const std::string& path);
private:
    PArUco::Params _params;
    PArUco::Detections _detections;
};

}
