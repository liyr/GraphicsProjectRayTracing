#pragma once
#include "Ray.h"

class Photon :
    public Ray
{
public:

    Photon(MyVector o_, MyVector d_): Ray(o_, d_)
    {
    }

    ~Photon()
    {
    }
};

