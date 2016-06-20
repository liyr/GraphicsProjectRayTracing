#pragma once
#include "Vector3d.h"

class Ray
{
public:

    MyVector RayPos, dir;

    Ray(MyVector o_, MyVector d_) : RayPos(o_), dir(d_)
    {
    }
};

