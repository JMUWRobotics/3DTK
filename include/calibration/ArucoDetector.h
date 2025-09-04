#ifndef CALIBRATION_ARUCODETECTOR_H
#define CALIBRATION_ARUCODETECTOR_H

#include "opencv2/core/version.hpp"

#if CV_MAJOR_VERSION > 3

#include "calibration/Detector.h"
#include "calibration/AprilTag.h"
#include <opencv2/opencv.hpp>
#include <opencv2/aruco.hpp>

#if (CV_MAJOR_VERSION >= 4) && (CV_MINOR_VERSION >= 7)
#include <opencv2/objdetect/aruco_detector.hpp>
#endif

namespace calibration
{

class ArucoDetector : public Detector
{
public:
    ArucoDetector(std::vector<AprilTag::AprilTag3f> patternPoints = std::vector<AprilTag::AprilTag3f>(), std::string dictionaryName = "DICT_6X6_250");
    ~ArucoDetector();

public:
    bool detect(const cv::Mat& image);
    void writeDetectionsToFile(const std::string& path);
    void readDetectionsFromFile(const std::string& path);

    static int dictionaryFromString(const std::string &str) noexcept(false);

private:
    std::vector<AprilTag::AprilTag3f> _patternPoints;
    std::string _dictionaryName;

#if (CV_MAJOR_VERSION >= 4) && (CV_MINOR_VERSION >= 7)
    cv::aruco::Dictionary _dictionary;
    cv::aruco::DetectorParameters _detectorParams;
    std::shared_ptr<cv::aruco::ArucoDetector> _detector;
#else
    cv::Ptr<cv::aruco::Dictionary> _dictionary;
    cv::Ptr<cv::aruco::DetectorParameters> _detectorParams;
#endif
    std::vector<AprilTag::AprilTag2f> _tags;
};

} // namespace calibration

#endif

#endif // CALIBRATION_ARUCODETECTOR_H
