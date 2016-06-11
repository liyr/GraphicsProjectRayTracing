#include "kdtree.h"
#include <algorithm>
#include <iostream>
#include <tuple>
#include <queue>
int dim = 0;

int cmp(const photon &a, const photon &b)
{
    if(a.pos[dim] > b.pos[dim])
    {
        return 0;
    }
    if (a.pos[dim] <= b.pos[dim])
    {
        return 1;
    }
    //if (a.pos[dim] < b.pos[dim])
    //{
    //    return -1;
    //}
}

kdtree::kdtree(std::vector<photon>* phoar, int rot, int ed, int dimen)
{
    phoarr = phoar;
    root = rot;
    end = ed;
    dim = dimen;
    if (dim == 3) dim = 0;
    std::sort(phoarr->begin() + root, phoarr->begin() + end, cmp);
    if (end - root <= 2) return;
    auto tmp = (*phoarr)[root];
    (*phoarr)[root] = (*phoarr)[root + (end - root)/2];
    (*phoarr)[root + (end - root) / 2] = tmp;
    std::sort(phoarr->begin() + root + 1, phoarr->begin() + end, cmp);
    index = root + (end - root) / 2 + 1;
    lchild = new kdtree(phoarr, root + 1, index - 1, dimen + 1);
    rchild = new kdtree(phoarr, index, ed, dimen + 1);
}

uchar* kdtree::knn(int i, int j, double cx, double cy, double cz)
{
    uchar* ans = new uchar[3];
    //std::priority_queue<photon,> kneighbor;
    ans[0] = 0;
    ans[1] = 0;
    ans[2] = 0;


    return ans;
}

kdtree::~kdtree()
{
    delete lchild, rchild;
}
