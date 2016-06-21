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
            switch(dim)
            {
            case 0:
                return abs((a.RayPos - pos).x) > abs((b.RayPos - pos).x);
            case 1:
                return abs((a.RayPos - pos).y) > abs((b.RayPos - pos).y);
            case 2:
                return abs((a.RayPos - pos).z) > abs((b.RayPos - pos).z);
            }
        }
    };
    MyVector knn(const MyVector pos_)
    {
        pos = pos_;
        std::priority_queue< Photon, std::vector<Photon>, cmp> pq;

        while(end - root >= 2)
        {
            switch(dimen)
            {
            case 0:
                break;
            case 1:
                break;
            case 2:
                break;
            }
            if((*phoarr)[root].RayPos.x > pos_.x)
            {
                
            }
        }


        MyVector res;
        const int k = 10;
        for (int i = 0; i < k; i++)
        {
            res = res + pq.top().color;
            pq.pop();
        }
        res = res * (1);
        return res;
    }
    ~kdtree()
    {
        delete lchild, rchild;
    }
};

