#include "calibration/CCTagDetector.h"

#include "cctag/Detection.hpp"
#include "opencv2/core/persistence.hpp"

namespace calibration {

using namespace cctag;

CCTagDetector::CCTagDetector(size_t nCrowns)
    : _params(nCrowns), _bank(nCrowns) {
    _objectPoints.clear();
}

bool CCTagDetector::detect(const cv::Mat& image_) {
    cv::Mat gray;
    CCTag::List markers;

    if (image_.channels() != 1)
        cvtColor(image_, gray, cv::COLOR_BGR2GRAY);
    else
        gray = image_;

    auto tick = std::chrono::high_resolution_clock::now();
    {
        cctagDetection(markers, 0, 0, gray, _params, _bank);
    }
    auto tock = std::chrono::high_resolution_clock::now();

    _detectionTime = std::chrono::duration_cast<std::chrono::milliseconds>(tock - tick).count();

    _imagePoints.clear();
    for (const auto &marker : markers)
        _imagePoints.push_back({ marker.x(), marker.y() });

    return _imagePoints.size() > 0;
}

void CCTagDetector::writeDetectionsToFile(const std::string& path) {
    cv::FileStorage file(path, cv::FileStorage::WRITE);
    const cv::Mat_<cv::Point2f> header(_imagePoints, false);

    file.write("imagePoints", header);
}

void CCTagDetector::readDetectionsFromFile(const std::string& path) {
    const cv::FileStorage file(path, cv::FileStorage::READ);
    cv::Mat_<cv::Point2f> data;

    file["imagePoints"] >> data;

    _imagePoints.clear();
    for (const auto &p : data)
        _imagePoints.push_back(p);
}

}