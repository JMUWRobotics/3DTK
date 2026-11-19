/*
 * sc_fixed nearest neighbor implementation
 *
 * Copyright (C) Fabia Frank
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief Implementation of fixed-point nearest neighbor search with kd-tree and range search
 * @author Fabia Frank. Institute of Computer Science, University of Wuerzburg, Germany.
 */



#include <vector>

#include "sc_fixed/sc_ICP.h"



struct bounding_box{
    fixed_val xmin, xmax, ymin, ymax, zmin, zmax;

    bounding_box(const std::vector<std::array<fixed_val, 3>>& points, int start, int end){
    if(end >= start){
    xmin = xmax = points[start][0];
    ymin = ymax = points[start][1];
    zmin = zmax = points[start][2];

    for(int i = start; i <= end; i++){
        if(points[i][0] < xmin) xmin = points[i][0];
        if(points[i][0] > xmax) xmax = points[i][0];
        if(points[i][1] < ymin) ymin = points[i][1];
        if(points[i][1] > ymax) ymax = points[i][1];
        if(points[i][2] < zmin) zmin = points[i][2];
        if(points[i][2] > zmax) zmax = points[i][2];
    }} else {
        xmin = xmax = ymin = ymax = zmin = zmax = 0;
    };    }
    bounding_box(){
        xmin = xmax = ymin = ymax = zmin = zmax = 0;

    };

};

struct box_leaf{
    bounding_box box;
    bool isLeaf;
    box_leaf(bounding_box box, bool isLeaf){
        this->box = box;
        this->isLeaf = isLeaf;
    };
};


   
struct kd_tree{
    std::vector<std::array<fixed_val, 3>> leafPoints;
    kd_tree* leftChild = nullptr;
    kd_tree* rightChild = nullptr;
    bool isLeaf= false;
    bounding_box box = bounding_box();
    int splitaxis;
    fixed_val splitvalue;

    public: kd_tree(std::vector<std::array<fixed_val, 3>>& points,int start, int end){
        if(end >= start){
            box = bounding_box(points, start, end);

            //stop-criterion amount of points
            if(end-start < 10){
                isLeaf = true;
                for(int i = start; i <= end; i++){
                    leafPoints.push_back(points[i]);
                }
                return;
            }

            //axis-length
            fixed_val dx = box.xmax - box.xmin;
            fixed_val dy = box.ymax - box.ymin;
            fixed_val dz = box.zmax - box.zmin;
            //stop-criterion distance of points
            if(dx < 0.01 && dy < 0.01 && dz < 0.01){
                isLeaf = true;
                for(int i = start; i <= end; i++){
                    leafPoints.push_back(points[i]);
                }
                return;
            }

            //Choose splitaxis (heuristic: longest extent of points along axis)
            if(dx > dy && dx > dz){
                splitaxis = 0;
                splitvalue = (box.xmin + box.xmax)/2;
            } else if (dy > dz){
                splitaxis = 1;
                splitvalue = (box.ymin + box.ymax)/2;
            }
            else{
                splitaxis = 2;
                splitvalue = (box.zmin + box.zmax)/2;
            }

        
            //partition of points
            int left = start;
            int right = end;
            while(left <= right){
                while (left <= end && points[left][splitaxis] < splitvalue){
                    left ++;
                }
                while (right >= start && points[right][splitaxis] >= splitvalue){
                right--;}
                if(left <= right){
                std::swap(points[left], points[right]);     //Prüfen ob umsortieren ok!!!!!!!!!!!!!!!!!!!!!
                left++;
                right--;
            }
            }
                leftChild = new kd_tree(points, start, right);
                rightChild = new kd_tree(points, right+1, end); 
        }

    }

    ~kd_tree(){
        delete leftChild;
        delete rightChild;}

    //Range search    
    void getPointsInRange(fixed_val range, std::array<fixed_val, 3> srcPoint, std::vector<std::array<fixed_val, 3>>&  pointsInRange){
        if( (box.xmax >= srcPoint[0] - range) && (box.xmin <= srcPoint[0] + range) &&
            (box.ymax >= srcPoint[1] - range) && (box.ymin <= srcPoint[1] + range) &&
            (box.zmax >= srcPoint[2] - range) && (box.zmin <= srcPoint[2] + range)) {
            if(isLeaf){
                for(unsigned int i = 0; i < leafPoints.size(); i++){
                        pointsInRange.push_back(leafPoints[i]);
                }
                return;
            } else{
                if(leftChild != nullptr) leftChild->getPointsInRange(range, srcPoint, pointsInRange);
                if(rightChild != nullptr) rightChild->getPointsInRange(range, srcPoint, pointsInRange);
            }}
        }
    

    bool nns(fixed_val mdm, std::vector<std::array<fixed_val, 3>>& source, std::vector<std::array<fixed_val, 3>>& matchedTarget, std::vector<std::array<fixed_val, 3>>& matchedSource){
        fixed_val sqr_mdm = mdm * mdm;
        for(unsigned int i = 0; i < source.size(); i++){
            std::array<fixed_val, 3> srcPoint = source[i];
            std::vector<std::array<fixed_val, 3>> pointsInRange;

            getPointsInRange(mdm, srcPoint, pointsInRange);

            if(pointsInRange.empty()){
                continue;
            }
            std::array<fixed_val, 3> min = pointsInRange[0];
            std::array<fixed_val, 3> tgtPoint;
            fixed_val sqr_dist;
            fixed_val sqr_dist_min = sqr_mdm;
            
                for(unsigned int j = 0; j < pointsInRange.size(); j++){
                    tgtPoint = pointsInRange[j];
                    fixed_val dx = srcPoint[0] - tgtPoint[0];
                    fixed_val dy = srcPoint[1] - tgtPoint[1];
                    fixed_val dz = srcPoint[2] - tgtPoint[2];
                    sqr_dist = dx * dx + dy * dy + dz * dz;
                    if(sqr_dist < sqr_mdm && sqr_dist < sqr_dist_min){
                    min = tgtPoint;
                    sqr_dist_min = sqr_dist;}
                }  
            
            if(sqr_dist_min < sqr_mdm){
                matchedSource.push_back(srcPoint);
                matchedTarget.push_back(min);
            }

    }
        if(matchedTarget.empty()) return false;
            return true;


   }
};
   bool nearest_neighbor(fixed_val mdm, std::vector<std::array<fixed_val, 3>>& source, std::vector<std::array<fixed_val, 3>>& target, std::vector<std::array<fixed_val, 3>>& matchedTarget, std::vector<std::array<fixed_val, 3>>& matchedSource){
        matchedSource.clear();
        matchedTarget.clear();
        if(source.empty() || target.empty()) return false;
        kd_tree tree = kd_tree(source,0, source.size() -1);
        return tree.nns(mdm, target, matchedSource, matchedTarget);
        
   }

