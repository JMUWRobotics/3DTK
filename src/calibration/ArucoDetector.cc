#include "calibration/ArucoDetector.h"
#include <chrono>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>
#include <calibration/DetectionFileHandler.h>

#if (CV_MAJOR_VERSION > 3)

namespace calibration {

int ArucoDetector::dictionaryFromString(const std::string &str) noexcept(false) {
    int dictionary_type;
    if (str.compare("DICT_6X6_250") == 0) {
        dictionary_type = cv::aruco::DICT_6X6_250;
    } else if (str.compare("DICT_6X6_1000") == 0) {
        dictionary_type = cv::aruco::DICT_6X6_1000;
    } else if (str.compare("DICT_APRILTAG_36h11") == 0) {
        dictionary_type = cv::aruco::DICT_APRILTAG_36h11;
    } else {
        throw(std::invalid_argument("unrecognized dictionary name, only DICT_6X6_250, DICT_6X6_1000 and DICT_APRILTAG_36h11 are supported!"));
    }
    return dictionary_type;
}

ArucoDetector::ArucoDetector(std::vector<AprilTag::AprilTag3f> patternPoints, std::string dictionaryName) :
    _patternPoints(patternPoints),
    _dictionaryName(dictionaryName)
{
    int dictionary_type = dictionaryFromString(dictionaryName);

    _dictionary = cv::aruco::getPredefinedDictionary(dictionary_type);
#if (CV_MAJOR_VERSION >= 4) && (CV_MINOR_VERSION >= 7)
    _detectorParams.cornerRefinementMethod = cv::aruco::CORNER_REFINE_SUBPIX;
    _detector = std::make_shared<cv::aruco::ArucoDetector>(_dictionary,_detectorParams);
#else
    _detectorParams = cv::aruco::DetectorParameters::create();
    _detectorParams->cornerRefinementMethod = cv::aruco::CORNER_REFINE_SUBPIX;
#endif

}

ArucoDetector::~ArucoDetector()
{
}

bool ArucoDetector::detect(const cv::Mat& image)
{
    _tags.clear();

    cv::Mat gray;

    if (image.channels() != 1) {
        cvtColor(image, gray, cv::COLOR_BGR2GRAY);
    } else {
        gray = image;
    }

    std::vector<int> markerIds;
    std::vector<std::vector<cv::Point2f> > markerCorners;

#if (CV_MAJOR_VERSION >= 4) && (CV_MINOR_VERSION >= 7)
    _detector->detectMarkers(gray, markerCorners, markerIds);
#else
    cv::aruco::detectMarkers(gray, _dictionary, markerCorners, markerIds, _detectorParams);
#endif

    for (size_t i = 0; i < markerIds.size(); i++) {
        AprilTag::AprilTag2f aprilTag2f = AprilTag::AprilTag2f(markerIds.at(i));

        aprilTag2f.point1 = markerCorners.at(i).at(0);
        aprilTag2f.point2 = markerCorners.at(i).at(1);
        aprilTag2f.point3 = markerCorners.at(i).at(2);
        aprilTag2f.point4 = markerCorners.at(i).at(3);

        _tags.push_back(aprilTag2f);
    }

    _imagePoints.clear();
    _objectPoints.clear();

    if (_patternPoints.size() == 0) {
        // if no pattern points provided copy all tags to image points
        for (AprilTag::AprilTag2f aprilTag2f : _tags) {
            _imagePoints.push_back(aprilTag2f.point1);
            _imagePoints.push_back(aprilTag2f.point2);
            _imagePoints.push_back(aprilTag2f.point3);
            _imagePoints.push_back(aprilTag2f.point4);
        }
    } else {
        // if pattern points provided copy only matching image points
        for (AprilTag::AprilTag2f aprilTag2f : _tags) {
            for (AprilTag::AprilTag3f aprilTag3f : _patternPoints) {
                if (aprilTag3f.compair(aprilTag2f) == 0) {
                    _imagePoints.push_back(aprilTag2f.point1);
                    _imagePoints.push_back(aprilTag2f.point2);
                    _imagePoints.push_back(aprilTag2f.point3);
                    _imagePoints.push_back(aprilTag2f.point4);

                    _objectPoints.push_back(aprilTag3f.point1);
                    _objectPoints.push_back(aprilTag3f.point2);
                    _objectPoints.push_back(aprilTag3f.point3);
                    _objectPoints.push_back(aprilTag3f.point4);
                }
            }
        }
    }

    return (_tags.size() > 0);
}

void ArucoDetector::writeDetectionsToFile(const std::string& path)
{
    std::map<std::string, cv::Point3f> detections;
    for(AprilTag::AprilTag2f tag : _tags){
        //corner 1
        std::stringstream id;
        id << "APR" << tag.id << "TL";
        detections.insert({id.str(), cv::Point3f(tag.point1.x, tag.point1.y, 0.0)});
        //corner 2
        id.str("");
        id << "APR" << tag.id << "BL";
        detections.insert({id.str(), cv::Point3f(tag.point2.x, tag.point2.y, 0.0)});
        //corner 3
        id.str("");
        id << "APR" << tag.id << "BR";
        detections.insert({id.str(), cv::Point3f(tag.point3.x, tag.point3.y, 0.0)});
        //corner 4
        id.str("");
        id << "APR" << tag.id << "TR";
        detections.insert({id.str(), cv::Point3f(tag.point4.x, tag.point4.y, 0.0)});
    }
    DetectionFileHandler::writeDetectionsToFile(path, detections);
}

void ArucoDetector::readDetectionsFromFile(const std::string& path)
{
    //TODO not compatible with the new detection file format
    std::cerr << "depricated function, don't use it! it's not compatible with the new detection file format!" << std::endl;
    _tags.clear();

    if (boost::filesystem::exists(boost::filesystem::path(path))) {
        _tags = AprilTag::createAprilTag2fFromFile(path);
    }
}

} // namespace calibration

#endif
