/*
 * feature implementation
 *
 * Copyright (C) Hamidreza Houshiar
 *
 * Released under the GPL version 3.
 *
 */

#include "slam6d/fbr/feature.h"

using namespace std;
using namespace cv;

namespace fbr {

    feature::feature() {
#ifdef WITH_OPENCV_NONFREE
        fDetectorMethod = SIFT_DET;
        fDescriptorMethod = SIFT_DES;
#else
        fDetectorMethod = ORB_DET;
        fDescriptorMethod = ORB_DES;
        fDetectorMethod = KAZE_DET;
        fDescriptorMethod = KAZE_DES;
        fDetectorMethod = AKAZE_DET;
        fDescriptorMethod = AKAZE_DES;
#endif
        fFiltrationMethod = DISABLE_FILTER;
    }

    feature::feature(feature_detector_method detector, feature_descriptor_method descriptor,
                     feature_filtration_method filtration) {
        fDetectorMethod = detector;
        fDescriptorMethod = descriptor;
        fFiltrationMethod = filtration;
    }

    void feature::featureDetection(cv::Mat pImage, feature_detector_method method, cv::Mat rImage,
                                   feature_filtration_method fMethod) {
        fDetectorMethod = method;
        fFiltrationMethod = fMethod;


        switch (fDetectorMethod) {
            //Detect the keypoints using SURF Detector
#ifdef WITH_OPENCV_NONFREE
            case SURF_DET:{
              double minHessian = 800;
#if (CV_MAJOR_VERSION >= 3) && (CV_MINOR_VERSION >= 0)
              Ptr<xfeatures2d::SURF> detector = xfeatures2d::SURF::create(minHessian);
              detector->detectAndCompute(pImage, noArray(), keypoints, descriptors);
#else
              cv::SurfFeatureDetector detector(minHessian);
              detector.detect(pImage, keypoints);
#endif
              break;
            }

              //Detect the keypoints using SIFT Detector
            case SIFT_DET:{
#if (CV_MAJOR_VERSION >= 4) && (CV_MINOR_VERSION >= 10)
              Ptr<cv::SIFT> detector = cv::SIFT::create();
              detector->detect(pImage, keypoints, descriptors);
#else
#if (CV_MAJOR_VERSION >= 3) && (CV_MINOR_VERSION >= 0)
              Ptr<xfeatures2d::SIFT> detector = xfeatures2d::SIFT::create();
              detector->detect(pImage, keypoints, descriptors);
#else
              cv::SiftFeatureDetector detector;
              detector.detect(pImage, keypoints);
#endif
#endif
              break;
            }
#endif
            //Detect the keypoints using ORB Detector
            case ORB_DET: {
#if CV_MAJOR_VERSION <= 2
                cv::OrbFeatureDetector detector;
                detector.detect(pImage, keypoints);
#else
                auto detector = cv::ORB::create();
                detector->detect(pImage, keypoints);
#endif
                break;
            }

                //Detect the keypoints using FAST Detector
            case FAST_DET: {
#if CV_MAJOR_VERSION <= 2
                cv::FastFeatureDetector detector;
                detector.detect(pImage, keypoints);
#else
                auto detector = cv::FastFeatureDetector::create();
                detector->detect(pImage, keypoints);
#endif
                break;
            }

                //Detect the keypoints using STAR Detector
#if CV_MAJOR_VERSION <= 2
            case STAR_DET: {
                cv::StarFeatureDetector detector;
                detector.detect(pImage, keypoints);

                break;
            }
#endif

                //Detect the keypoints using KAZE Detector
#if (CV_MAJOR_VERSION >= 3) && (CV_MINOR_VERSION >= 0)
            case KAZE_DET: {

                Ptr<cv::KAZE> detector = cv::KAZE::create();
                detector->detect(pImage, keypoints);

                break;
            }
#endif

                //Detect the keypoints using AKAZE Detector
#if (CV_MAJOR_VERSION >= 3) && (CV_MINOR_VERSION >= 0)
            case AKAZE_DET: {

                Ptr<cv::AKAZE> detector = cv::AKAZE::create();
                detector->detect(pImage, keypoints);
            }
#endif
        }

        featureFiltration(pImage, rImage);
    }


    void feature::featureDetection(cv::Mat pImage, feature_detector_method method) {
        cv::Mat rImage;
        featureDetection(pImage, method, rImage, fFiltrationMethod);
    }

    void feature::featureDetection(cv::Mat pImage) {
        featureDetection(pImage, fDetectorMethod);
    }

    void feature::featureDescription(cv::Mat pImage, feature_descriptor_method method) {

        fDescriptorMethod = method;
        if (keypoints.size() == 0)
            featureDetection(pImage);
        switch (fDescriptorMethod) {
#ifdef WITH_OPENCV_NONFREE
            case SURF_DES:{
              //Create descriptor using SURF

#if (CV_MAJOR_VERSION >= 3) && (CV_MINOR_VERSION >= 0)
              Ptr<xfeatures2d::SURF> extractor = xfeatures2d::SURF::create();
              extractor->compute(pImage, keypoints, descriptors);
#else
              cv::SurfDescriptorExtractor extractor;
              extractor.compute(pImage, keypoints, descriptors);
#endif
              break;
            }
            case SIFT_DES:{
              //Create descriptor using SIFT
#if (CV_MAJOR_VERSION >= 4) && (CV_MINOR_VERSION >= 10)
              Ptr<cv::SIFT> detector = cv::SIFT::create();
              detector->compute(pImage, keypoints, descriptors);
#else
#if (CV_MAJOR_VERSION >= 3) && (CV_MINOR_VERSION >= 0)
              Ptr<xfeatures2d::SIFT> extractor = xfeatures2d::SIFT::create();
              extractor->compute(pImage, keypoints, descriptors);
#else
              cv::SiftDescriptorExtractor extractor;
              extractor.compute(pImage, keypoints, descriptors);
#endif
#endif
              break;
            }
#endif
            case ORB_DES: {
                //Create descriptor using ORB
#if CV_MAJOR_VERSION <= 2
                cv::OrbDescriptorExtractor extractor;
                extractor.compute(pImage, keypoints, descriptors);
#else
                auto extractor = cv::ORB::create();
                extractor->compute(pImage, keypoints, descriptors);
#endif
                break;
            }
            case KAZE_DES: {
#if CV_MAJOR_VERSION >= 3 && CV_MINOR_VERSION >= 0
                Ptr<cv::KAZE> detector = cv::KAZE::create();
                detector->compute(pImage, keypoints, descriptors);
                break;
#endif
            }
            case AKAZE_DES: {
#if CV_MAJOR_VERSION >= 3 && CV_MINOR_VERSION >= 0
                Ptr<cv::AKAZE> detector = cv::AKAZE::create();
                detector->compute(pImage, keypoints, descriptors);
                break;
#endif
            }
        }
    }

    void feature::featureDescription(cv::Mat pImage) {
        featureDescription(pImage, fDescriptorMethod);
    }

    feature_detector_method feature::getDetectorMethod() {
        return fDetectorMethod;
    }

    feature_descriptor_method feature::getDescriptorMethod() {
        return fDescriptorMethod;
    }

    feature_filtration_method feature::getFeatureFiltrationMethod() {
        return fFiltrationMethod;
    }

    //check for the keypoints vector not to be empty
    vector <cv::KeyPoint> feature::getFeatures() {
        return keypoints;
    }

    //check for the descriptor Mat not to be empty
    cv::Mat feature::getDescriptors() {
        return descriptors;
    }

    void feature::setFeatures(vector <cv::KeyPoint> keypoint) {
        keypoints = keypoint;
    }

    void feature::setDescriptors(cv::Mat descriptor) {
        descriptors = descriptor;
    }

    void feature::getDescription(description_method method) {
        if (method == FEATURE_DESCRIPTION)
            cout << "fDetectorMethod: " << featureDetectorMethodToString(fDetectorMethod)
                 << ", number of detected features: " << keypoints.size() << ", feature filtration method: "
                 << featureFiltrationMethodToString(fFiltrationMethod) << "." << endl;
        else if (method == DESCRIPTOR_DESCRIPTION)
            cout << "fDescriptorMethod: " << featureDescriptorMethodToString(fDescriptorMethod) << "." << endl;
        else
            cout << "fDetectorMethod: " << featureDetectorMethodToString(fDetectorMethod)
                 << ", number of detected features: " << keypoints.size() << ", feature filtration method: "
                 << featureFiltrationMethodToString(fFiltrationMethod) << ", fDescriptorMethod: "
                 << featureDescriptorMethodToString(fDescriptorMethod) << "." << endl;
        cout << endl;
    }

    unsigned int feature::getNumberOfFeatures() {
        return keypoints.size();
    }

    void feature::featureFiltration(cv::Mat pImage, cv::Mat rImage) {
        vector <cv::KeyPoint> filteredKeypoints;
        if (fFiltrationMethod == OCCLUSION) {
            for (unsigned int i = 0; i < keypoints.size(); i++) {
                int x, y;
                x = keypoints[i].pt.x;
                y = keypoints[i].pt.y;
                float range[8];
                if (rImage.at<float>(y, x) != 0) {
                    range[0] = rImage.at<float>(y, x) - rImage.at<float>(y + 1, x + 1);
                    range[1] = rImage.at<float>(y, x) - rImage.at<float>(y, x + 1);
                    range[2] = rImage.at<float>(y, x) - rImage.at<float>(y - 1, x + 1);
                    range[3] = rImage.at<float>(y, x) - rImage.at<float>(y - 1, x);
                    range[4] = rImage.at<float>(y, x) - rImage.at<float>(y - 1, x - 1);
                    range[5] = rImage.at<float>(y, x) - rImage.at<float>(y, x - 1);
                    range[6] = rImage.at<float>(y, x) - rImage.at<float>(y + 1, x - 1);
                    range[7] = rImage.at<float>(y, x) - rImage.at<float>(y + 1, x);

                    int count = 0;
                    for (unsigned int j = 0; j < 8; j++) {
                        if (range[j] < 20)
                            count++;
                    }
                    if (count == 8)
                        filteredKeypoints.push_back(keypoints[i]);
                }
            }
        } else if (fFiltrationMethod == STANDARD_DEVIATION) {
            for (unsigned int i = 0; i < keypoints.size(); i++) {
                int x, y;
                x = keypoints[i].pt.x;
                y = keypoints[i].pt.y;
                float range[9];
                if (rImage.at<float>(y, x) != 0) {
                    range[0] = rImage.at<float>(y, x) - rImage.at<float>(y + 1, x + 1);
                    range[1] = rImage.at<float>(y, x) - rImage.at<float>(y, x + 1);
                    range[2] = rImage.at<float>(y, x) - rImage.at<float>(y - 1, x + 1);
                    range[3] = rImage.at<float>(y, x) - rImage.at<float>(y - 1, x);
                    range[4] = rImage.at<float>(y, x) - rImage.at<float>(y - 1, x - 1);
                    range[5] = rImage.at<float>(y, x) - rImage.at<float>(y, x - 1);
                    range[6] = rImage.at<float>(y, x) - rImage.at<float>(y + 1, x - 1);
                    range[7] = rImage.at<float>(y, x) - rImage.at<float>(y + 1, x);
                    range[8] = rImage.at<float>(y, x) - rImage.at<float>(y, x);

                    double t_std, r = 0, temp = 0;
                    for (unsigned int j = 0; j < 9; j++)
                        r += range[j];
                    r /= 9;
                    for (unsigned int j = 0; j < 9; j++)
                        temp += ((range[j] - r) * (range[j] - r));
                    t_std = sqrt(temp / 9);
                    if (t_std < 0.1)
                        filteredKeypoints.push_back(keypoints[i]);
                }
            }
        }
        if (fFiltrationMethod != DISABLE_FILTER)
            keypoints = filteredKeypoints;
    }
}

