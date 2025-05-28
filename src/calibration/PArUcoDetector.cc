#include "calibration/PArUcoDetector.h"

#include "calibration/ArucoDetector.h"

namespace calibration {

PArUcoDetector::PArUcoDetector(float tagExpandScale, const std::string& dictionaryName)
    : _params()
{
    _objectPoints.clear();

    _params.aruco.dictionary = cv::aruco::getPredefinedDictionary(ArucoDetector::dictionaryFromString(dictionaryName));
    _params.tagExpandScale = tagExpandScale;
}

bool PArUcoDetector::detect(const cv::Mat& image_) {
    cv::Mat gray;

    _imagePoints.clear();

    if (image_.channels() != 1)
        cvtColor(image_, gray, cv::COLOR_BGR2GRAY);
    else
        gray = image_;

    auto tick = std::chrono::high_resolution_clock::now();
    {
        PArUco::detect(gray, _detections, _params);
    }
    auto tock = std::chrono::high_resolution_clock::now();

    _detectionTime = std::chrono::duration_cast<std::chrono::milliseconds>(tock - tick).count();

    for (const auto &d : _detections) {
        for (const auto &p : d.circleCenters) {
            if (!p.has_value()) continue;

            _imagePoints.push_back(*p);
        }
    }

    return _imagePoints.size() > 0;
}

void PArUcoDetector::writeDetectionsToFile(const std::string& path) {
    cv::FileStorage file(path, cv::FileStorage::WRITE);
    const cv::Mat_<cv::Point2f> header(_imagePoints, false);

    file.write("imagePoints", header);
}

void PArUcoDetector::readDetectionsFromFile(const std::string& path) {
    const cv::FileStorage file(path, cv::FileStorage::READ);
    cv::Mat_<cv::Point2f> data;

    file["imagePoints"] >> data;

    _imagePoints.clear();
    for (const auto &p : data)
        _imagePoints.push_back(p);
}

}