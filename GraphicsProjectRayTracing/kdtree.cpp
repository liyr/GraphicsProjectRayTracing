#include "kdtree.h"
#include <algorithm>
#include "operation.h"
#include <iostream>
#include <tuple>
#include <queue>
#include "proces.h"
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
    double tmp[3] = { 0,0,0 };
    double pos[3];
    pos[0] = -500 + i;
    pos[1] = 500 - j;
    pos[2] = 500;
    normalize(pos);
    int index = -1;
    double dis;
    photon trace;
    double look[3] = { 500,500,0 };
    setValue(trace.pos, look);
    setValue(trace.dir, pos);
    findIntersec(trace, index, dis);
    trace.pos[0] += trace.dir[0] * dis;
    trace.pos[1] += trace.dir[1] * dis;
    trace.pos[2] += trace.dir[2] * dis;
    double minusres[3];// , crossres[3];
    for(auto it = phoarr->begin(); it != phoarr->end(); ++it)
    {
        setValue(minusres, it->pos);
        minus(minusres, trace.pos);
        //cross_product(minusres, pos, crossres);
        auto x = norm(minusres);
        tmp[0] += 0.1 * exp(-x*x/2000)*it->p[0];
        tmp[1] += 0.1 * exp(-x*x/2000)*it->p[1];
        tmp[2] += 0.1 * exp(-x*x/2000)*it->p[2];
    }
    ans[0] = static_cast<int>(tmp[0]);
    ans[1] = static_cast<int>(tmp[1]);
    ans[2] = static_cast<int>(tmp[2]);
    return ans;
}

kdtree::~kdtree()
{
    delete lchild, rchild;
}
