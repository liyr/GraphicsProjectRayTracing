#pragma once
#include "Object.h"
#include <vector>

class Scene
{
public:
    std::vector<Object*> objects;

    Scene()
    {
    }

    bool intersect(const Ray& r, double& t, int& id)
    {
        double d;
        t = INFINITY;
        for (auto it = objects.begin(); it != objects.end(); ++it)
            if ((d = (*it)->intersect(r)) && d < t)
            {
                t = d;
                id = it - objects.begin();
            }
        return t < INFINITY;
    }


    ~Scene()
    {
    }
};

