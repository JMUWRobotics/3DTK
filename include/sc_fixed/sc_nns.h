#ifndef SC_NNS_H
#define SC_NNS_H

#include <array>
#include <vector>
#include "sc_fixed/sc_ICP.h"


struct bounding_box{
    fixed_val xmin, xmax, ymin, ymax, zmin, zmax;

    bounding_box(const std::vector<std::array<fixed_val, 3>>& points, int start, int end);
    bounding_box();};


struct box_leaf{
    bounding_box box;
    bool isLeaf;
    box_leaf(bounding_box box, bool isLeaf);
};


struct kd_tree{
    std::vector<std::array<fixed_val, 3>> leafPoints;
    kd_tree* leftChild = nullptr;
    kd_tree* rightChild = nullptr;
    bool isLeaf= false;
    bounding_box box = bounding_box();
    int splitaxis;
    fixed_val splitvalue;

    kd_tree(std::vector<std::array<fixed_val, 3>>& points,int start, int end);

    ~kd_tree();

    void getPointsInRange(fixed_val range, std::array<fixed_val, 3> srcPoint, std::vector<std::array<fixed_val, 3>>&  pointsInRange);
    bool nns(fixed_val mdm, std::vector<std::array<fixed_val, 3>>& source, std::vector<std::array<fixed_val, 3>>& matchedTarget, std::vector<std::array<fixed_val, 3>>& matchedSource);
    };


bool nearest_neighbor(fixed_val mdm, std::vector<std::array<fixed_val, 3>>& source, std::vector<std::array<fixed_val, 3>>& target, std::vector<std::array<fixed_val, 3>>& matchedTarget, std::vector<std::array<fixed_val, 3>>& matchedSource);

#endif