#pragma once
#include <core/cvdef.h>

class Space
{
public:
    int depth, length, width;

    Space(int ,int ,int);
    ~Space();
};

class photon
{
public:
    double pos[3] = {0,0,0};       // position ( 3 x 32 bit floats )
    double p[3]; 
    double dir[3];
    double phi = 0, theta = 0;    
    short flag = 0;          // flag used for kd-tree
    static photon emitPhoton(int = 500,int = 500,int = 0, int = 1, int = -1, int = -1);
    static photon backLight();
};

class Camera
{
public:
    int x, y, z;

};

