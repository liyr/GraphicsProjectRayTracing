#pragma once
#include "Ray.h"

class Photon :
    public Ray
{
public:
    MyVector color;
    Photon(MyVector o_, MyVector d_, MyVector c_): Ray(o_, d_), color(c_)
    {
    }

    ~Photon()
    {
    }
};

