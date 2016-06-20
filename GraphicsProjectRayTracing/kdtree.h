#pragma once
#include <vector>
#include <queue>
#include "Photon.h"
#include <algorithm>
#include <core/cvdef.h>
int dim = 0;
MyVector pos;
int cmpa(const Photon &a, const Photon &b)
{
    switch(dim)
    {
    case 0:
        return a.RayPos.x > b.RayPos.x ? 0 : 1;
        break;
    case 1:
        return a.RayPos.y > b.RayPos.y ? 0 : 1;
        break;
    case 2:
        return a.RayPos.z > b.RayPos.z ? 0 : 1;
        break;
    }
}



class kdtree
{
public:
    int root, end, dimen, index;
    kdtree *lchild, *rchild;
    std::vector<Photon>* phoarr;
    kdtree(std::vector<Photon>* phoar, int rot, int ed, int dimen)
    {
        phoarr = phoar;
        root = rot;
        end = ed;
        dim = dimen;
        if (dim == 3) dim = 0;
        std::sort(phoarr->begin() + root, phoarr->begin() + end, cmpa);
        if (end - root <= 2) return;
        auto tmp = (*phoarr)[root];
        (*phoarr)[root] = (*phoarr)[root + (end - root) / 2];
        (*phoarr)[root + (end - root) / 2] = tmp;
        std::sort(phoarr->begin() + root + 1, phoarr->begin() + end, cmpa);
        index = root + (end - root) / 2 + 1;
        lchild = new kdtree(phoarr, root + 1, index - 1, dimen + 1);
        rchild = new kdtree(phoarr, index, ed, dimen + 1);
    }
    struct cmp {
        bool operator()(Photon a, Photon b) const
        {
            return (a.RayPos - pos).norm() > (b.RayPos - pos).norm();
        }
    };
    MyVector knn(const MyVector pos_)
    {
        pos = pos_;
        std::priority_queue< Photon, std::vector<Photon>, cmp> pq;



        MyVector res;
        return res;
    }
    ~kdtree()
    {
        delete lchild, rchild;
    }
};

