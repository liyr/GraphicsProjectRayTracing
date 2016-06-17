#include "Space.h"
#include <random>
#include <cmath>
#include "operation.h"
#include <iostream>

std::uniform_real_distribution<double> disz(0, 1);
std::uniform_real_distribution<double> disx(-1, 1);
std::uniform_real_distribution<double> disy(-1, 1);
std::default_random_engine generator;

Space::Space(int length, int width, int depth): depth(depth), length(length), width(width)
{

}

Space::~Space()
{
}

photon photon::emitPhoton(int cx, int cy, int cz, int num, int aimx, int aimy)
{
    photon ans;
    
    double light[3] = { cx, cy, cz };
    if (aimx == -1 || aimy == -1) {
        do
        {
            ans.pos[0] = disx(generator);
            ans.pos[1] = disy(generator);
            ans.pos[2] = disz(generator);
        } while (norm(ans.pos) > 1);
    }
    else
    {
        ans.pos[0] = -500 + aimx;
        ans.pos[1] = 500 - aimy; 
        ans.pos[2] = 500;
    }
    //double tmp[3];

    normalize(ans.pos);
    //setValue(tmp, ans.pos);
    setValue(ans.dir, ans.pos);
    //normalize(ans.dir);

    //if(ans.dir[1] != ans.pos[1])
    //{
    //    std::cout << "djfshgfk" << std::endl;
    //}
    ans.theta = acos(ans.pos[2]);
    ans.phi = atan2(ans.pos[0], ans.pos[1]);
    plus(ans.pos, light);
    ans.p[0] = 255. / num / 4;
    ans.p[1] = 255. / num / 4;
    ans.p[2] = 255. / num / 4;
    return ans;
}

photon photon::backLight()
{
    photon ans;
    //if (rand() % 2 == 0) {
    //    ans.pos[0] = disx(generator) * 500;
    //    ans.pos[1] = disy(generator) * 500;
    //    ans.pos[2] = 0;
    //    ans.p[0] = 255;
    //    ans.p[1] = 255;
    //    ans.p[2] = 255;
    //    randomHalfSphere(ans.theta, ans.phi, ans.dir);
    //    ans.dir[0] = sin(ans.theta) * sin(ans.phi);
    //    ans.dir[1] = sin(ans.theta) * cos(ans.phi);
    //    ans.dir[2] = cos(ans.theta);
    //}
    //else
    //{
        ans.p[0] = 255;
        ans.p[1] = 255;
        ans.p[2] = 255;
        randomHalfSphere(ans.theta, ans.phi, ans.dir);
        double centr[3];
        centr[0] = 300;
        centr[1] = 1020;
        centr[2] = 1000;

        ans.dir[0] = sin(ans.theta) * sin(ans.phi);
        ans.dir[1] = -cos(ans.theta);
        ans.dir[2] = sin(ans.theta) * cos(ans.phi);

        ans.pos[0] = centr[0] + 100 * ans.dir[0];
        ans.pos[1] = centr[0] + 100 * ans.dir[1];
        ans.pos[2] = centr[0] + 100 * ans.dir[2];

    //}
    return ans;
}
