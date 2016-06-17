#pragma once
#include <iostream>  
#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>
#include "Space.h"
#include "kdtree.h"
#include "operation.h"
#include <random>
#include <ctime>
#include <thread>
using namespace cv;
using point = double[3];
static const int length = 1000, width = 1000;


static const int photonMapIteTimes = 50000, objectNum = 12;
static const int depth = 1000;
static const int cx = 500, cy = 500, cz = 0, lz = 500;
static point zero = { 0,0,0 };
static const double fluxThreshold = 1. / 1000;
static const int ittimes = 5120;

static double *brdf_chrome, *brdf_delrin, *brdf_hematite, *brdf_natural, *brdf_walnut, *brdf_telfon, *brdf_violet;
static std::vector<photon> phoarray;

static Mat bg_mat;

static Mat marble_mat;

static Mat wd_mat;

static Mat sk_mat;

static Mat eh_mat;

inline void findIntersec(const photon &, int &, double &);


inline bool sphereIntersec(const photon& pho, double& ans, double pos2centr[3], double centr[3], double radi, double& value1)
{
    setValue(pos2centr, centr);
    minus(pos2centr, pho.pos);
    auto l = norm(pos2centr);
    auto dot = dot_product(pos2centr, pho.dir);
    //if (dot < 0) return false;
    if (radi * radi - l*l + dot * dot >= -1e-4)
    {
        ans = ((dot - sqrt(radi * radi - l*l + dot * dot + 1e-4)) > 1e-5) ?
            (dot - sqrt(radi * radi - l*l + dot * dot + 1e-4)) :
            (dot + sqrt(radi * radi - l*l + dot * dot + 1e-4));
        if (ans < 1e-5) return false;
        value1 = ans;
        return true;
    }
    return false;
}

inline bool checkValue(const photon& pho, double ans, double& value1)
{
    if (ans > 1 && pho.pos[0] + ans * pho.dir[0] >= -1e-2 && pho.pos[0] + ans * pho.dir[0] <= 1000 + 1e-2
        && pho.pos[1] + ans * pho.dir[1] >= -1e-2 && pho.pos[1] + ans * pho.dir[1] <= 1000 + 1e-2
        && pho.pos[2] + ans * pho.dir[2] >= 500 - 1e-2 && pho.pos[2] + ans * pho.dir[2] <= 1500 + 1e-2)
    {
        value1 = ans;
        return true;
    }
    return false;
}

inline double intersec(const photon &pho, int ind)// , int iup, int idown)
{
    double ans;

    double pos2centr[3];
    double centr[3] = { 500,100,1400 };
    double radi = 200;
    double value1;


    switch (ind)
    {
    case 0:
        ans = (1000 + lz - pho.pos[2]) / pho.dir[2];
        if (checkValue(pho, ans, value1)) return value1;
        break;
    case 1:
        ans = (1000 - pho.pos[1]) / pho.dir[1];
        if (checkValue(pho, ans, value1)) return value1;
        break;
    case 2:
        ans = (-pho.pos[1]) / pho.dir[1];
        if (checkValue(pho, ans, value1)) return value1;
        break;
    case 3:
        ans = (1000 - pho.pos[0]) / pho.dir[0];
        if (checkValue(pho, ans, value1)) return value1;
        break;
    case 4:
        //if(pho.pos[0] - pho.dir[0] != 500 || pho.pos[1] - pho.dir[1] != 500 || pho.pos[2] != pho.dir[2])
        //{
        //    std::cout << "poaurihr" << std::endl;
        //}
        ans = (-pho.pos[0]) / pho.dir[0];
        if (checkValue(pho, ans, value1)) return value1;
        //if (ans > 0) return ans;
        break;
    case 5:
        centr[0] = 400;
        centr[1] = 200;
        centr[2] = 1250;
        radi = 100;
        if (sphereIntersec(pho, ans, pos2centr, centr, radi, value1)) return value1;
        break;
    case 6:
        centr[0] = 300;
        centr[1] = 250;
        centr[2] = 1000;
        radi = 100;
        if (sphereIntersec(pho, ans, pos2centr, centr, radi, value1)) return value1;
        break;
    case 7:
        centr[0] = 700;
        centr[1] = 200;
        centr[2] = 950;
        radi = 200;
        if (sphereIntersec(pho, ans, pos2centr, centr, radi, value1)) return value1;
        break;
    case 8:
        centr[0] = 100;
        centr[1] = 250;
        centr[2] = 750;
        radi = 100;
        if (sphereIntersec(pho, ans, pos2centr, centr, radi, value1)) return value1;
        break;
    case 9:
        //centr[0] = 900;
        //centr[1] = 250;
        //centr[2] = 700;
        //radi = 100;
        //if (sphereIntersec(pho, ans, pos2centr, centr, radi, value1)) return value1;
        break;
    case 10:
        ans = (-pho.pos[2]) / pho.dir[2];
        if (ans > 1) return ans;
        break;
    case 11:
        centr[0] = 300;
        centr[1] = 1020;
        centr[2] = 1000;
        radi = 100;
        if (sphereIntersec(pho, ans, pos2centr, centr, radi, value1)) return value1;
        break;
    }

    return INFINITY;
}

inline void calBiNormal(double* normal, double* bi_normal)
{
    if (abs(normal[0]) > 1e-2)
    {
        bi_normal[0] = -normal[1];
        bi_normal[1] = normal[0];
        bi_normal[2] = 0;
    }
    else
    {
        bi_normal[0] = 0;
        bi_normal[1] = -normal[2];
        bi_normal[2] = normal[1];
    }
}

inline void calNewNormal(double normal[3], double bi_normal[3], int index, double pos[3])
{
    double centr[3];
    switch (index)
    {
    case 0:
        normal[0] = 0, normal[1] = 0, normal[2] = -1;
        bi_normal[0] = 0, bi_normal[1] = 1, bi_normal[2] = 0;
        break;
    case 1:
        normal[0] = 0, normal[1] = -1, normal[2] = 0;
        bi_normal[0] = 1, bi_normal[1] = 0, bi_normal[2] = 0;
        break;
    case 2:
        normal[0] = 0, normal[1] = 1, normal[2] = 0;
        bi_normal[0] = 1, bi_normal[1] = 0, bi_normal[2] = 0;
        break;
    case 3:
        normal[0] = -1, normal[1] = 0, normal[2] = 0;
        bi_normal[0] = 0, bi_normal[1] = 1, bi_normal[2] = 0;
        break;
    case 4:
        normal[0] = 1, normal[1] = 0, normal[2] = 0;
        bi_normal[0] = 0, bi_normal[1] = 1, bi_normal[2] = 0;
        break;
    case 5:
        centr[0] = 400;
        centr[1] = 200;
        centr[2] = 1250;
        setValue(normal, pos);
        minus(normal, centr);
        normalize(normal);
        calBiNormal(normal, bi_normal);
        break;
    case 6:
        centr[0] = 300;
        centr[1] = 250;
        centr[2] = 1000;
        setValue(normal, pos);
        minus(normal, centr);
        normalize(normal);
        calBiNormal(normal, bi_normal);
        break;
    case 7:
        centr[0] = 700;
        centr[1] = 200;
        centr[2] = 950;
        setValue(normal, pos);
        minus(normal, centr);
        normalize(normal);
        calBiNormal(normal, bi_normal);
        break;
    case 8:
        centr[0] = 100;
        centr[1] = 250;
        centr[2] = 750;
        setValue(normal, pos);
        minus(normal, centr);
        normalize(normal);
        calBiNormal(normal, bi_normal);
        break;
    case 9:
        centr[0] = 900;
        centr[1] = 250;
        centr[2] = 700;
        setValue(normal, pos);
        minus(normal, centr);
        normalize(normal);
        calBiNormal(normal, bi_normal);
        break;
    case 10:
        normal[0] = 0, normal[1] = 0, normal[2] = 1;
        bi_normal[0] = 1, bi_normal[1] = 0, bi_normal[2] = 0;
        break;
    case 11:
        centr[0] = 300;
        centr[1] = 1020;
        centr[2] = 1000;
        setValue(normal, pos);
        minus(normal, centr);
        normalize(normal);

        break;
    }
}

inline void moveTo(photon &pho, int index, double &ans)
{
    point tmp;
    setValue(tmp, pho.dir);
    double normal[3], bi_normal[3];
    pho.pos[0] += ans * tmp[0];
    pho.pos[1] += ans * tmp[1];
    pho.pos[2] += ans * tmp[2];
    calNewNormal(normal, bi_normal, index, pho.pos);
    refresh(pho, normal, bi_normal);
    uchar* data;
    int atmp;
    switch (index)
    {
    case 0:
        //ans = (1000 + lz - pho.pos[2]) / tmp[2];
        //pho.p[0] *= 0.5;
        data = marble_mat.ptr<uchar>(int(pho.pos[1]));
        pho.p[0] *= (int)data[3 * (int)pho.pos[0]] / 255.;
        pho.p[1] *= (int)data[3 * (int)pho.pos[0] + 1] / 255.;
        pho.p[2] *= (int)data[3 * (int)pho.pos[0] + 2] / 255.;
        break;
    case 1:
        //ans = (1000 - pho.pos[1]) / tmp[1];
        //pho.p[1] *= 0.15;
        data = sk_mat.ptr<uchar>(int(pho.pos[0]));
        pho.p[0] *= (int)data[3 * (int)(pho.pos[2] - 500)] / 255.;
        pho.p[1] *= (int)data[3 * (int)(pho.pos[2] - 500) + 1] / 255.;
        pho.p[2] *= (int)data[3 * (int)(pho.pos[2] - 500) + 2] / 255.;
        break;
    case 2:
        //ans = - pho.pos[1] / tmp[1];
        //pho.p[2] *= 0.15;
        data = wd_mat.ptr<uchar>(int(pho.pos[0]));
        pho.p[0] *= (int)data[3 * (int)(pho.pos[2] - 500)] / 255.;
        pho.p[1] *= (int)data[3 * (int)(pho.pos[2] - 500) + 1] / 255.;
        pho.p[2] *= (int)data[3 * (int)(pho.pos[2] - 500) + 2] / 255.;
        break;
    case 3:
        pho.p[0] *= 0.75;
        //ans = (1000 - pho.pos[0]) / tmp[0];
        //data = bg_mat.ptr<uchar>(int(pho.pos[1]));
        //pho.p[0] *= (int)data[3 * (int)(pho.pos[2] - 500)] / 255.;
        //pho.p[1] *= (int)data[3 * (int)(pho.pos[2] - 500) + 1]/ 255.;
        //pho.p[2] *= (int)data[3 * (int)(pho.pos[2] - 500) + 2]/ 255.;
        break;
    case 4:
        pho.p[0] *= 0.25;
        pho.p[1] *= 0.25;
        //ans = - pho.pos[0] / tmp[0];
        //data = bg_mat.ptr<uchar>(int(pho.pos[1]));
        //pho.p[0] *= (int)data[3 * (int)(pho.pos[2] - 500)]/ 255.;
        //pho.p[1] *= (int)data[3 * (int)(pho.pos[2] - 500) + 1]/ 255.;
        //pho.p[2] *= (int)data[3 * (int)(pho.pos[2] - 500) + 2]/ 255.;
        break;
    case 5:
        break;
    case 6:
        //pho.p[0] *= 0.15;
        //pho.p[1] *= 0.85;
        break;
    case 7:
        //pho.p[0] *= 0.15;
        //pho.p[1] *= 0.85;
        data = eh_mat.ptr<uchar>(int(1024 - pho.pos[1] * 2.56));
        atmp = (int)(atan2(pho.pos[2] - 950, pho.pos[0] - 700) * 1024 / M_PI + 1024);
        //std::cout << atmp << std::endl;
        //if(atmp < 0 || atmp > 2047)
        //std::cout << "fuck" << std::endl;
        pho.p[0] *= ((int)data[3 * atmp]) / 255.;
        pho.p[1] *= ((int)data[3 * atmp + 1]) / 255.;
        pho.p[2] *= ((int)data[3 * atmp + 2]) / 255.;
        break;
    case 8:
        pho.p[1] *= 0.85;
        break;
    case 9:
        pho.p[1] *= 0.15;
        break;
    case 10:
        //ans = -pho.pos[2] / tmp[2];
        pho.p[0] *= 0.50;// .2;//abs(cos(pho.theta));
        pho.p[1] *= 0.50;//.2;//abs(cos(pho.theta));
        pho.p[2] *= 0.50;// .2;//abs(cos(pho.theta));
        break;
    case 11:
        //ans = -pho.pos[2] / tmp[2];
        pho.p[0] *= 60;// abs(cos(pho.theta));
        pho.p[1] *= 60;// abs(cos(pho.theta));
        pho.p[2] *= 60;// abs(cos(pho.theta));
        break;
    }
    pho.p[0] *= 0.9;//abs(cos(pho.theta));
    pho.p[1] *= 0.9;//abs(cos(pho.theta));
    pho.p[2] *= 0.9;//abs(cos(pho.theta));
}

inline int Russia(photon &pho, int index)
{
    //if (rand() % 100 == 0) return -1;
    return 1;
}


inline photon genNewPho(photon &pho, int index)
{
    double normal[3], bi_normal[3], tri_normal[3];
    calNewNormal(normal, bi_normal, index, pho.pos);
    cross_product(normal, bi_normal, tri_normal);
    double theta;
    double phi;
    point pos;
    switch (index)
    {
    case 0:
        break;
    case 2:
        break;
    case 3:
        break;
    case 4:
        break;
    case 5:
        theta = M_PI - pho.theta;
        phi = pho.phi;
        break;
    case 6:
        if (pho.theta > M_PI / 2) {
            if (sin(pho.theta) > 1 / 1.44) theta = M_PI - pho.theta;
            else theta = M_PI - asin(sin(pho.theta) * 1.44);
        }
        else
        {
            if (sin(pho.theta) > 1 / 1.44) theta = pho.theta;
            else theta = asin(sin(pho.theta) * 1.44);
        }
        phi = pho.phi;
        break;
    case 8:
        break;
    default:
        randomHalfSphere(theta, phi, pos);
        break;
    }
    photon ans = pho;

    ans.theta = theta;
    ans.phi = phi;
    ans.dir[0] = cos(ans.theta) * normal[0] + sin(ans.theta) * cos(ans.phi) * bi_normal[0] + sin(ans.theta) * sin(ans.phi) * tri_normal[0];
    ans.dir[1] = cos(ans.theta) * normal[1] + sin(ans.theta) * cos(ans.phi) * bi_normal[1] + sin(ans.theta) * sin(ans.phi) * tri_normal[1];
    ans.dir[2] = cos(ans.theta) * normal[2] + sin(ans.theta) * cos(ans.phi) * bi_normal[2] + sin(ans.theta) * sin(ans.phi) * tri_normal[2];
    //if (index == 5)
    //{
    //    std::cout << pho.dir[0] <<" "<< pho.dir[1] << " " << pho.dir[2] << std::endl;
    //    std::cout << ans.dir[0] << " " << ans.dir[1] << " " << ans.dir[2] << std::endl;
    //    std::cout << ans.theta << " " << ans.phi << std::endl;
    //    std::cout << "fuck" << std::endl;
    //}
    return ans;
}

inline double strength(photon &trace)
{
    return trace.p[0] + trace.p[1] + trace.p[2];
}

inline void LightTracing()
{
    photon pho = photon::backLight();
    while (1)
    {
        if (strength(pho) < fluxThreshold) { break; }
        int index = -1;
        double ans;
        findIntersec(pho, index, ans);
        if (index == -1) break;
        moveTo(pho, index, ans);
        if (index != 5 && index != 6) phoarray.push_back(pho);
        if (Russia(pho, index) == -1) break;
        if (index == 10 || index == 11) break;
        pho = genNewPho(pho, index);
    }
}


inline void findIntersec(const photon &trace, int& index, double &ans)
{
    double dep = INFINITY;
    for (int it = 0; it < objectNum; ++it)
    {
        double curDep = intersec(trace, it);
        if (curDep < dep)
        {
            dep = curDep;
            index = it;
        }
    }
    ans = dep;
}

inline void rayT(int it, int jt, point tmp, int i)
{
    auto trace = photon::emitPhoton(500, 500, 0, i, it, jt);
    while (1)
    {
        if (strength(trace) < fluxThreshold) { setValue(trace.p, zero); break; }
        int index = -1;
        double ans;
        findIntersec(trace, index, ans);
        if (index == -1)
        {
            setValue(trace.p, zero);
            break;
        }
        //if (index == 10 || index == 11 || Russia(trace, index) == -1) break;
        moveTo(trace, index, ans);
        if (index == 10 || index == 11 || Russia(trace, index) == -1) break;
        //if (index == 11) break;
        //break;
        trace = genNewPho(trace, index);
    }
    plus(tmp, trace.p);
}

inline uchar* rayTracing(int it, int jt, int cx, int cy, int cz)
{
    uchar* ans = new uchar[3];
    point tmp = { 0,0,0 };

    for (int i = 0; i < ittimes; ++i) {
        rayT(it, jt, tmp, ittimes / 4);
    }
    ans[0] = min(static_cast<int>(tmp[0]), 255);
    ans[1] = min(static_cast<int>(tmp[1]), 255);
    ans[2] = min(static_cast<int>(tmp[2]), 255);
    return ans;
}

inline void colorPaint(Mat& input, kdtree &root)
{
    int rowNumber = input.rows;
    int colNumber = input.cols * input.channels();

    for (int i = 0; i < rowNumber; ++i)
    {
        if (i % 4 == 0)std::cout << "now cal the row " << i << std::endl;
        uchar* data = input.ptr<uchar>(i);
        uchar* al = 0;
        uchar* ar = 0;
        for (int j = 0; j < colNumber; ++j)
        {
            if (j % 3 == 0)
            {
                //delete al; al = root.knn(j / 3, i, 500, 500, 0);
                delete ar; ar = rayTracing(j/3, i, 500, 500, 0); 
            }
            data[j] = min((int)ar[j % 3], 255);
        }
    }
}

inline void print(double vec[3])
{
    std::cout << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
}


inline void preCalBrdf()
{
    read_brdf("./model/chrome.binary", brdf_chrome);
    read_brdf("./model/delrin.binary", brdf_delrin);
    read_brdf("./model/hematite.binary", brdf_hematite);
    read_brdf("./model/natural-209.binary", brdf_natural);
    read_brdf("./model/special-walnut-224.binary", brdf_walnut);
    read_brdf("./model/teflon.binary", brdf_telfon);
    read_brdf("./model/violet-acrylic.binary", brdf_violet);
}

