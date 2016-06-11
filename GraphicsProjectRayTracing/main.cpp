#include <iostream>  
#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>
#include "Space.h"
#include "kdtree.h"
#include "operation.h"
#include <random>

using namespace cv;
const int times = 10, num = 11;
const int depth = 1000;
const int length = 1000;
const int width = 1000;
const int cx = 500, cy = 500, cz = 0, lz = 500;
double zero[3] = { 0,0,0 };

std::vector<photon> phoarray;
std::uniform_real_distribution<double> mdisz(0, 1);
std::uniform_real_distribution<double> mdisx(-1, 1);
std::uniform_real_distribution<double> mdisy(-1, 1);
std::default_random_engine mgenerator;


double intersec(const photon &pho, int ind)
{
    double tmp[3];
    tmp[0] = pho.dir[0];
    tmp[1] = pho.dir[1];
    tmp[2] = pho.dir[2];
    double ans;
    switch(ind)
    {
    case 0:
        ans = (1000 + lz - pho.pos[2]) / tmp[2];
        if (ans > 0) return ans;
        break;
    case 1:
        ans = (1000 - pho.pos[1]) / tmp[1];
        if (ans > 0) return ans;
        break;
    case 2:
        ans = - pho.pos[1] / tmp[1];
        if (ans > 0) return ans;
        break;
    case 3:
        ans = (1000 - pho.pos[0]) / tmp[0];
        if (ans > 0) return ans;
        break;
    case 4:
        ans = - pho.pos[0] / tmp[0];
        if (ans > 0) return ans;
        break;
    case 5:
        break;
    case 6:
        break;
    case 7:
        break;
    case 8:
        break;
    case 9:
        break;
    case 10:
        ans = - pho.pos[2] / tmp[2];
        if (ans > 0) return ans;
        break;
    }

    return INFINITY;
}

void calNewNormal(double normal[3], double bi_normal[3], int index, double pos[3])
{
    switch(index)
    {
    case 0:
        normal[0] = 0, normal[1] = 0, normal[2] = -1;
        bi_normal[0] = 0, bi_normal[1] = 1, bi_normal[2] = 0;
        break;
    case 1:
        normal[0] = -1, normal[1] = 0, normal[2] = 0;
        bi_normal[0] = 0, bi_normal[1] = 1, bi_normal[2] = 0;
        break;
    case 2:
        normal[0] = 1, normal[1] = 0, normal[2] = 0;
        bi_normal[0] = 0, bi_normal[1] = 1, bi_normal[2] = 0;
        break;
    case 3:
        normal[0] = 0, normal[1] = -1, normal[2] = 0;
        bi_normal[0] = 1, bi_normal[1] = 0, bi_normal[2] = 0;
        break;
    case 4:
        normal[0] = 0, normal[1] = 1, normal[2] = 0;
        bi_normal[0] = 1, bi_normal[1] = 0, bi_normal[2] = 0;
        break;
    case 5:
        break;
    case 6:
        break;
    case 7:
        break;
    case 8:
        break;
    case 9:
        break;
    case 10:
        break;
    case 11:
        normal[0] = 0, normal[1] = 0, normal[2] = 1;
        bi_normal[0] = 1, bi_normal[1] = 0, bi_normal[2] = 0;
        break;
    }
}

void moveTo(photon &pho, int index)
{
    double tmp[3];
    tmp[0] = pho.dir[0];
    tmp[1] = pho.dir[1];
    tmp[2] = pho.dir[2];
    double ans;
    double normal[3], bi_normal[3];
    switch(index)
    {
    case 0:
        ans = (1000 + lz - pho.pos[2]) / tmp[2];
        pho.pos[0] += ans * tmp[0];
        pho.pos[1] += ans * tmp[1];
        pho.pos[2] += ans * tmp[2];
        calNewNormal(normal, bi_normal, index, pho.pos);
        refresh(pho, normal, bi_normal);
        pho.p[0] *= sin(pho.phi) * 0.5;
        pho.p[1] *= sin(pho.phi);
        pho.p[2] *= sin(pho.phi);
        break;
    case 1:
        ans = (1000 - pho.pos[1]) / tmp[1];
        pho.pos[0] += ans * tmp[0];
        pho.pos[1] += ans * tmp[1];
        pho.pos[2] += ans * tmp[2];
        calNewNormal(normal, bi_normal, index, pho.pos);
        refresh(pho, normal, bi_normal);
        pho.p[0] *= sin(pho.phi);
        pho.p[1] *= sin(pho.phi) * 0.5;
        pho.p[2] *= sin(pho.phi);
        break;
    case 2:
        ans = - pho.pos[1] / tmp[1];
        pho.pos[0] += ans * tmp[0];
        pho.pos[1] += ans * tmp[1];
        pho.pos[2] += ans * tmp[2];
        calNewNormal(normal, bi_normal, index, pho.pos);
        refresh(pho, normal, bi_normal);
        pho.p[0] *= sin(pho.phi);
        pho.p[1] *= sin(pho.phi);
        pho.p[2] *= sin(pho.phi) * 0.5;
        break;
    case 3:
        ans = (1000 - pho.pos[0]) / tmp[0];
        pho.pos[0] += ans * tmp[0];
        pho.pos[1] += ans * tmp[1];
        pho.pos[2] += ans * tmp[2];
        calNewNormal(normal, bi_normal, index, pho.pos);
        refresh(pho, normal, bi_normal);
        pho.p[0] *= sin(pho.phi) * 0.5;
        pho.p[1] *= sin(pho.phi) * 0.5;
        pho.p[2] *= sin(pho.phi);
        break;
    case 4:
        ans = - pho.pos[0] / tmp[0];
        pho.pos[0] += ans * tmp[0];
        pho.pos[1] += ans * tmp[1];
        pho.pos[2] += ans * tmp[2];
        calNewNormal(normal, bi_normal, index, pho.pos);
        refresh(pho, normal, bi_normal);
        pho.p[0] *= sin(pho.phi) * 0.5;
        pho.p[1] *= sin(pho.phi);
        pho.p[2] *= sin(pho.phi) * 0.5;
        break;
    case 5:
        break;
    case 6:
        break;
    case 7:
        break;
    case 8:
        break;
    case 9:
        break;
    case 10:
        // need done
        ans = -pho.pos[2] / tmp[2];
        pho.pos[0] += ans * tmp[0];
        pho.pos[1] += ans * tmp[1];
        pho.pos[2] += ans * tmp[2];
        calNewNormal(normal, bi_normal, index, pho.pos);
        refresh(pho, normal, bi_normal);
        pho.p[0] *= sin(pho.phi);
        pho.p[1] *= sin(pho.phi);
        pho.p[2] *= sin(pho.phi);
        break;
    }
}

int Russia(photon &pho, int index)
{
    //int cho = -1;
    if (rand() % 100 == 0) return -1;
    return 1;
}

photon genNewPho(photon pho, int index)
{
    double normal[3], bi_normal[3], tri_normal[3];
    calNewNormal(normal, bi_normal, index, pho.pos);
    cross_product(normal, bi_normal, tri_normal);

    double pos[3];
    do
    {
        pos[0] = mdisx(mgenerator);
        pos[1] = mdisy(mgenerator);
        pos[2] = mdisz(mgenerator);
    } while (norm(pos) > 1);
    normalize(pos);
    auto theta = acos(pos[2]);
    auto phi = atan2(pos[0], pos[1]);

    photon ans = pho;



    ans.theta = theta;
    ans.phi = phi;
    ans.dir[0] += cos(ans.theta) * normal[0] + sin(ans.theta) * cos(ans.phi) * bi_normal[0] + sin(ans.theta) * sin(ans.phi) * tri_normal[0];
    ans.dir[1] += cos(ans.theta) * normal[1] + sin(ans.theta) * cos(ans.phi) * bi_normal[1] + sin(ans.theta) * sin(ans.phi) * tri_normal[1];
    ans.dir[2] += cos(ans.theta) * normal[2] + sin(ans.theta) * cos(ans.phi) * bi_normal[2] + sin(ans.theta) * sin(ans.phi) * tri_normal[2];
    return ans;
}

void LightTracing(photon pho)
{
    while (1)
    {
        double depth = INFINITY;
        int index = -1;
        for (int i = 0; i < num; ++i)
        {
            double curDep = intersec(pho, i);
            if (curDep < depth)
            {
                depth = curDep;
                index = i;
            }
        }
        if (index == -1) break;
        moveTo(pho, index);
        phoarray.push_back(pho);
        if (Russia(pho, index) == -1) break;
        pho = genNewPho(pho,index);
    }
}

bool atLightSource(photon &trace)
{
    return trace.pos[2] < 1e-2;
}

double strength(photon &trace)
{
    return trace.p[0] + trace.p[1] + trace.p[2];
}

uchar* rayTracing(int it, int jt, int cx, int cy, int cz)
{
    const double HOLD = 1. / 1000;
    uchar* ans = new uchar[3];
    ans[0] = 0;
    ans[1] = 0;
    ans[2] = 0;
    double tmp[3] = { 0,0,0 };
    const int ittimes = 255;
    for (int i = 0; i < ittimes; ++i) {
        auto trace = photon::emitPhoton(500,500,0,255,it,jt);
        while (1)
        {
            double depth = INFINITY;
            int index = -1;
            for (int it = 0; it < num; ++it)
            {
                double curDep = intersec(trace, it);
                if (curDep < depth)
                {
                    depth = curDep;
                    index = it;
                }
            }
            if (index == -1 || strength(trace) < HOLD) { trace.p[0] = 0, trace.p[1] = 0, trace.p[2] = 0; break; }
            moveTo(trace, index);

            if (atLightSource(trace) || Russia(trace, index) == -1) break;
            trace = genNewPho(trace, index);
        }
        tmp[0] += trace.p[0];
        tmp[1] += trace.p[1];
        tmp[2] += trace.p[2];
    }
    ans[0] = static_cast<int>(tmp[0]);
    ans[1] = static_cast<int>(tmp[1]);
    ans[2] = static_cast<int>(tmp[2]);
    return ans;
}

void colorPaint(Mat& input, kdtree &root)
{
    int rowNumber = input.rows;
    int colNumber = input.cols * input.channels();

    for(int i = 0; i < rowNumber; ++i)
    {
        uchar* data = input.ptr<uchar>(i);
        uchar* al = 0;
        uchar* ar = 0;
        for (int j = 0; j < colNumber; ++j)
        {
            if (j % 3 == 0) { delete al; al = root.knn(i, j, 500, 500, 0); delete ar; ar = rayTracing(i, j/3, 500, 500, 0); }
            data[j] = al[j % 3] + ar[j % 3];
        }
    }
}

void print(double vec[3])
{
    std::cout << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
}


int main()
{
	Mat img(length, width, CV_8UC3, Scalar::all(255));
	randu(img, Scalar::all(0), Scalar::all(255));
    Space space(length, width, depth);

    double *brdf_chrome, *brdf_delrin, *brdf_hematite, *brdf_natural, *brdf_walnut, *brdf_telfon, *brdf_violet;
    read_brdf("./model/chrome.binary", brdf_chrome);
    read_brdf("./model/delrin.binary", brdf_delrin);
    read_brdf("./model/hematite.binary", brdf_hematite);
    read_brdf("./model/natural-209.binary", brdf_natural);
    read_brdf("./model/special-walnut-224.binary", brdf_walnut);
    read_brdf("./model/teflon.binary", brdf_telfon);
    read_brdf("./model/violet-acrylic.binary", brdf_violet);



    while(1)
    {
        for (int i = 0; i < times; ++i)
        {
            if (rand() % 2 == 0)
            {
                auto a = photon::emitPhoton();
            }
            else
            {
                photon a = photon::backLight();
            }
            //lighttracing(phoarray[i]);
        }
        //TODO:save the map
        break;
        
    }
    //for (int i = 0; i < times; ++i)
    //    std::cout << phoarray[i].pos[0] << " " << phoarray[i].pos[1] << " " << phoarray[i].pos[2] << std::endl;
	//TODO: add objects(five sphere? + 1 bunny?) + Creature Texture
	//TODO: background light
	//TODO: (ray-tracing + photon mapping), multi-texture, soft shadows,Depth of field，mirroring, Anti-Aliasing 
    kdtree root_kdtree(&phoarray, 0, phoarray.end() - phoarray.begin(), 0);
    
    colorPaint(img, root_kdtree);


   //test:
    //double a[3] = { 0,0,1 };
    //double b[3] = { 0,1,0 };
    //double c[3] = { 1,0,0 };
    //double d[3] = { rand(),rand(),rand()};
    //print(d);
    //normalize(d);
    //print(d);
    //double comoving[3];
    //cross_product(a, b, c);
    //print(c);
    //comoving[2] = dot_product(d, a);
    //comoving[1] = dot_product(d, b);
    //comoving[0] = mix_product(d, a, b);
    //auto theta = acos(comoving[2]);
    //auto phi = atan2(comoving[0], comoving[1]);
    //std::cout << theta * 180 / M_PI << " " << phi * 180 / M_PI << std::endl;
    //double ans[3] = { 0,0,0 };
    //print(ans);
    //std::cout << std::endl;
    //for (int i = 0; i < times; ++i)
    //    std::cout << phoarray[i].pos[0] << " " << phoarray[i].pos[1] << " " << phoarray[i].pos[2] << std::endl;
	namedWindow("效果示例");
	imshow("效果示例", img); 
	waitKey(6000);
	imwrite(".\\ans.png", img);

    system("pause");
	return 0;
}
