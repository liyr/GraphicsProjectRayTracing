#pragma once
#include "Object.h"
#include "Vector3d.h"
#include "Ray.h"

class Rectan :
    public Object
{
public:
    const MyVector normal;
    const MyVector e, color;

    Rectan(MyVector normal_, MyVector bi_normal_, MyVector p_, MyVector e_, MyVector c_, Refl_t refl_) :
        normal(normal_), e(e_), color(c_), Object(p_, refl_)
    {
    }

    double intersect(const Ray& r) const override
    {
        double dis = (pos - r.RayPos).dot(r.dir);
        if (dis < 0 || !(r.RayPos + r.dir * dis).inbox())return 0;
        return dis;
    }

    MyVector getColor(const MyVector& pos) const override
    {
        return color;
    }

    MyVector getNormal(const MyVector& pos) const override
    {
        return normal;
    }

    MyVector getEmission(const MyVector& pos) const override
    {
        return e;
    }

    ~Rectan()
    {
    }
};

