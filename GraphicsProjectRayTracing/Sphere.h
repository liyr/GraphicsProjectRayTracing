#pragma once
#include "Object.h"
#include "Vector3d.h"
#include "Ray.h"
#include <random>
#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>

std::normal_distribution<double> dist(0.00, 0.2);
std::default_random_engine gen;


class Sphere :
    public Object
{
public:
    double rad; // radius
    MyVector e, color; // position, emission, color
    cv::Mat wood_mat = cv::imread("./wooden-floor-texture-wallpaper-10.jpg");

    Sphere(double rad_, MyVector p_, MyVector e_, MyVector c_, Refl_t refl_) :
        rad(rad_), e(e_), color(c_), Object(p_, refl_)
    {
    }

    double intersect(const Ray& r) const override
    { // returns distance, 0 if nohit
        MyVector op = pos - r.RayPos;
        double t, eps = 1e-4, b = op|r.dir, det = b * b - (op|op) + rad * rad;
        if (det < 0) return 0;
        else det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }

    MyVector getColor(const MyVector& position) const override
    {
        if (refl == WOOD)
        {
            auto x = (position - pos) + MyVector(100, 0, 100);
            x.x = std::min(500., x.x);
            x.x = std::max(.0, x.x);
            x.z = std::min(500., x.z);
            x.z = std::max(.0, x.z);
            const uchar* data = wood_mat.ptr<uchar>(int(x.x * 2));
            //std::cout << x.x << std::endl;
            MyVector tmp;
            tmp.x = ((int)data[3 * (int)(x.z * 2)]) / 255.;
            tmp.y = ((int)data[3 * (int)(x.z * 2) + 1]) / 255.;
            tmp.z = ((int)data[3 * (int)(x.z * 2) + 2]) / 255.;
            return tmp;
        }
        else
        {
            return color;
        }
    }

    MyVector getNormal(const MyVector& x) const override
    {
        if (refl == BUMP || refl == WARD)
        {
            MyVector t1((MyVector(0, 1, 0) % (x - pos)) * 0.2 * sin(((x - pos).x + (x - pos).y + (x - pos).z)/2));
            return (x - pos + t1).normalize();
        }
        else if (refl == WOOD)
        {
            MyVector t1(dist(gen), dist(gen), dist(gen));
            return (x - pos + t1).normalize();
        }
        else
            return (x - pos).normalize();
    }

    MyVector getEmission(const MyVector& pos) const override
    {
        return e;
    }

};

