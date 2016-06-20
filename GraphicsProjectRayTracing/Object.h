#pragma once
#include "Vector3d.h"
class Ray;

enum Refl_t
{
    DIFF,
    SPEC,
    REFR,
    WARD,
    WOOD,
    BUMP,
    BRDF
};


class Object
{
public:
    Refl_t refl;
    MyVector pos;

    Object(MyVector pos_, Refl_t refl_): refl(refl_), pos(pos_)
    {
    }

    virtual double intersect(const Ray& r) const = 0;

    virtual MyVector getColor(const MyVector& pos) const = 0;

    virtual MyVector getNormal(const MyVector& pos) const = 0;

    virtual MyVector getEmission(const MyVector& pos) const = 0;

};

