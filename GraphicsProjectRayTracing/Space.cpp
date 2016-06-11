#include "Space.h"
#include <random>
#include <cmath>
#include "operation.h"

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
    if(aimx == -1 && aimy == -1)
    do
    {
        ans.pos[0] = disx(generator);
        ans.pos[1] = disy(generator);
        ans.pos[2] = disz(generator);
    } while (norm(ans.pos) > 1);
    else
    {
        ans.pos[0] = -500 + aimx;
        ans.pos[1] = 500 - aimy; 
        ans.pos[2] = 500;
    }


    normalize(ans.pos);
    ans.dir[0] = ans.pos[0];
    ans.dir[1] = ans.pos[1];
    ans.dir[2] = ans.pos[2];
    ans.theta = acos(ans.pos[2]);
    ans.phi = atan2(ans.pos[0], ans.pos[1]);
    plus(ans.pos, light);
    ans.p[0] = 255./num;
    ans.p[1] = 255./num;
    ans.p[2] = 255./num;
    return ans;
}

photon photon::backLight()
{
    photon ans;
    ans.pos[0] = disx(generator) * 500;
    ans.pos[1] = disy(generator) * 500;
    ans.pos[2] = 0;
    ans.theta = 0;
    ans.phi = 0;
    ans.p[0] = 255;
    ans.p[1] = 255;
    ans.p[2] = 255;
    ans.dir[0] = 0;
    ans.dir[1] = 0;
    ans.dir[2] = 1;
    return ans;
}
