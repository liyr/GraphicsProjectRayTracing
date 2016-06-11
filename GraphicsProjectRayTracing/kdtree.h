#pragma once
#include "Space.h"
#include <vector>
#include <core/cvdef.h>

class kdtree
{
public:
    int root, end, dimen, index;
    kdtree *lchild, *rchild;
    std::vector<photon>* phoarr;
    kdtree(std::vector<photon>* phoar, int, int, int);
    uchar* knn(int, int, double, double, double);
    ~kdtree();
};

