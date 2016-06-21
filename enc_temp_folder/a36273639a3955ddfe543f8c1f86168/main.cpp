#pragma comment(linker, "/STACK:1024000000,1024000000") 

#include <iostream> 
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp>
#include <random>
#include <ctime>
#include "Ray.h"
#include "Sphere.h"
#include <functional>
#include <cmath>
#include "operation.h"
#include "Photon.h"
#define M_PI	3.1415926535897932384626433832795

using namespace cv;

double truncate(double x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x)
{
    return int(pow(truncate(x), 1 / 2.2) * 255 + .5);
}

double* brdf_m;

bool intersect(const Ray& r, double& t, int& id, const std::vector<Sphere>& spheres)
{
    double d;
    t = INFINITY;
    for (auto it = spheres.begin(); it != spheres.end(); ++it)
        if ((d = it->intersect(r)) && d < t)
        {
            t = d;
            id = it - spheres.begin();
        }
    return t < INFINITY;
}

MyVector getBiNormal(MyVector w)
{
    auto u = ((fabs(w.x) > .1 ? MyVector(0, 1) : MyVector(1)) % w).normalize();
    return u;
}

void sphereDiffusion(MyVector nl, double r1, double r2, MyVector& d)
{
    double r2s = sqrt(r2);
    MyVector w = nl;
    MyVector u;
    MyVector v;
    u = getBiNormal(w);
    v = w % u;
    d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize();
}

void wardDiffusion(MyVector in, MyVector nl, double r1, double r2, MyVector& d)
{
    MyVector w = nl;
    MyVector u;
    MyVector v;
    u = getBiNormal(w);
    v = w % u;
    double alpha_y = 0.1, alpha_x = 0.1;
    double phi = atan(alpha_y / alpha_x * tan(2 * M_PI * r2));
    double theta = atan(sqrt(-log(r1) / (cos(phi) * cos(phi) / alpha_x / alpha_x
        + sin(phi) * sin(phi) / alpha_y / alpha_y)));
    d = w * cos(theta) + v * sin(theta) * cos(phi) + u * sin(theta) * sin(phi);
    d = in - d * 2 * (d | in);
}

void brdf(const Ray& r, MyVector nl, MyVector& f, MyVector d)
{
    double theta_out, phi_out, theta_in, phi_in;
    r.dir.calSphereRepre(nl, getBiNormal(nl), theta_out, phi_out);
    d.calSphereRepre(nl, getBiNormal(nl), theta_in, phi_in);
    MyVector tmp;
    lookup_brdf_val(brdf_m, theta_in, phi_in, M_PI - theta_out, phi_out, tmp.x, tmp.y, tmp.z);
    f = f * tmp * 20;
}

void phong(const Ray ray, MyVector nl, MyVector &color, MyVector d)
{
    double kd = 1, ks = 0.7;
    color = color * (kd * (ray.dir | nl) + ks * pow((d | ray.dir),512));
}

MyVector radiance(const Ray& r, int depth, const std::vector<Sphere>& spheres, const std::function<double (void)>& ran)
{
    double t; // distance to intersection
    int id = 0; // id of intersected object
    if (!intersect(r, t, id, spheres)) return MyVector(); // if miss, return black
    const Object& obj = spheres[id]; // the hit object
    MyVector x = r.RayPos + r.dir * t, n = obj.getNormal(x), nl = (n | (r.dir)) < 0 ? n : n * -1, f = obj.getColor(x);
    double p = max(max(f.x, f.y), f.z); // max refl
    if (++depth > 5) if (ran() < p) f = f * (1 / p); else return obj.getEmission(x); //R.R.
    double r1, r2, r3;
    MyVector d;
    switch (obj.refl)
    {
    case DIFF:
    case BUMP:
        r1 = 2 * M_PI * ran() , r2 = ran();
        sphereDiffusion(nl, r1, r2, d);
        return obj.getEmission(x) + f * (radiance(Ray(x, d), depth, spheres, ran));
        break;
    case BRDF:
        r1 = 2 * M_PI * ran() , r2 = ran();
        sphereDiffusion(nl, r1, r2, d);
        brdf(r, nl, f, d);
        return obj.getEmission(x) + f * (radiance(Ray(x, d), depth, spheres, ran));
        break;
    case PHONG:
        r1 = 2 * M_PI * ran(), r2 = ran();
        sphereDiffusion(nl, r1, r2, d);
        phong(r, nl, f, d);
        return obj.getEmission(x) + f * (radiance(Ray(x, d), depth, spheres, ran));
        break;
    case WARD:
        r1 = ran() , r2 = ran();
        wardDiffusion(r.dir, nl, r1, r2, d);
        return obj.getEmission(x) + f * (radiance(Ray(x, d), depth, spheres, ran));
        break;
    case WOOD:
        r1 = 2 * M_PI * ran() , r2 = ran();
        sphereDiffusion(nl, r1, r2, d);
        return obj.getEmission(x) + f * (ran() > 0.65 ?
               radiance(Ray(x, d), depth, spheres, ran) : 
            radiance(Ray(x, r.dir - n * 2 * (n | (r.dir))), depth, spheres, ran));
        break;
    case SPEC:
        return obj.getEmission(x) + f * (radiance(Ray(x, r.dir - n * 2 * (n | (r.dir))), depth, spheres, ran));
        break;
    case REFR:
        Ray reflRay(x, r.dir - n * 2 * (n | (r.dir))); // Ideal dielectric REFRACTION
        bool into = (n | nl) > 0; // Ray from outside going in?
        double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.dir | (nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
            return obj.getEmission(x) + f * (radiance(reflRay, depth, spheres, ran));


        MyVector tdir = (r.dir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
        double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir | (n));
        double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
        return obj.getEmission(x) + f * (depth > 2 ? (ran() < P ? // Russian roulette
               radiance(reflRay, depth, spheres, ran) * RP : radiance(Ray(x, tdir), depth, spheres, ran) * TP) :
               radiance(reflRay, depth, spheres, ran) * Re + radiance(Ray(x, tdir), depth, spheres, ran) * Tr);
        break;
    }
}

Photon randomPhoton(const std::vector<Sphere>& spheres, const std::function<double(void)>& ran)
{
    double x, z;
    do
    {
        x = ran();
        z = ran();
    }
    while ((x - 0.5) * (x - 0.5) + (z - 0.5) * (z - 0.5) > 0.25) ;
    double r1 = 2 * M_PI * ran(), r2 = ran();
    MyVector dir;
    sphereDiffusion(MyVector(0, -1, 0), r1, r2, dir);
    Photon res(MyVector(36 * (x - 0.5) + 50, 81.6, 36 * (z - 0.5) + 81.6), dir, MyVector(1., 1., 1.));
    return res;
}

void PhotonTrace(int depth, const std::vector<Sphere>& spheres, const std::function<double(void)>& ran)
{
    Photon r = randomPhoton(spheres, ran);
    double t;
    int id = 0; // id of intersected object
    if (!intersect(r, t, id, spheres)) return; // if miss, return black
    const Object& obj = spheres[id]; // the hit object
    MyVector x = r.RayPos + r.dir * t, n = obj.getNormal(x), nl = (n | (r.dir)) < 0 ? n : n * -1, f = obj.getColor(x);
    double p = max(max(f.x, f.y), f.z); // max refl
    if (++depth > 5) if (ran() < p) f = f * (1 / p); else return; //R.R.
    double r1, r2, r3;
    MyVector d;
    switch (obj.refl)
    {
    case DIFF:
    case BUMP:
        r1 = 2 * M_PI * ran() , r2 = ran();
        sphereDiffusion(nl, r1, r2, d);
        obj.getEmission(x) + f * (radiance(Ray(x, d), depth, spheres, ran));
        break;
    case BRDF:
        r1 = 2 * M_PI * ran() , r2 = ran();
        sphereDiffusion(nl, r1, r2, d);
        brdf(r, nl, f, d);
        obj.getEmission(x) + f * (radiance(Ray(x, d), depth, spheres, ran));
        break;
    case WARD:
        r1 = ran() , r2 = ran();
        wardDiffusion(r.dir, nl, r1, r2, d);
        obj.getEmission(x) + f * (radiance(Ray(x, d), depth, spheres, ran));
        break;
    case WOOD:
        r1 = 2 * M_PI * ran() , r2 = ran();
        sphereDiffusion(nl, r1, r2, d);
        obj.getEmission(x) + f * (ran() > 0.65 ?
                                             (radiance(Ray(x, d), depth, spheres, ran)) : 
            (radiance(Ray(x, r.dir - n * 2 * (n | (r.dir))), depth, spheres, ran)));
        break;
    case SPEC:
        obj.getEmission(x) + f * (radiance(Ray(x, r.dir - n * 2 * (n | (r.dir))), depth, spheres, ran));
        break;
    case REFR:
        Ray reflRay(x, r.dir - n * 2 * (n | (r.dir))); // Ideal dielectric REFRACTION
        bool into = (n | nl) > 0; // Ray from outside going in?
        double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.dir | (nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
            obj.getEmission(x) + f * (radiance(reflRay, depth, spheres, ran));


        MyVector tdir = (r.dir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
        double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir | (n));
        double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
        obj.getEmission(x) + f * (depth > 2 ? (ran() < P ? // Russian roulette
        radiance(reflRay, depth, spheres, ran) * RP : radiance(Ray(x, tdir), depth, spheres, ran) * TP) :
        radiance(reflRay, depth, spheres, ran) * Re + radiance(Ray(x, tdir), depth, spheres, ran) * Tr);
        break;
    }
}

void emitPhotons(const std::vector<Sphere>& spheres, const std::function<double(void)>& ran)
{
    const int num = 10;
    for (int i = 0; i < num; ++i)
    {
        PhotonTrace(0, spheres, ran);
    }
}


int main(int argc, char* argv[])
{
    auto begin = clock();
    int w = 1024, h = 768, samps = 2000; // # samples
    Ray cam(MyVector(50, 52, 295.6), MyVector(0, -0.042612, -1).normalize()); // cam pos, dir
    MyVector cx = MyVector(w * .5135 / h), cy = (cx % cam.dir).normalize() * .5135, r, *c = new MyVector[w * h];
    read_brdf("./model/chrome.binary", brdf_m);
    const std::vector<Sphere> spheres = {//Scene: radius, position, emission, color, material
        Sphere(1e5, MyVector(1e5 + 1, 40.8, 81.6), MyVector(), MyVector(.95, .05, .05), BRDF),//Left
        Sphere(1e5, MyVector(-1e5 + 99, 40.8, 81.6), MyVector(), MyVector(.05, .05, .95), WARD),//Rght
        Sphere(1e5, MyVector(50, 40.8, 1e5), MyVector(), MyVector(.05, .05, .85), PHONG),//Back
        Sphere(1e5, MyVector(50, 40.8, -1e5 + 170), MyVector(), MyVector(), DIFF),//Frnt
        Sphere(1e5, MyVector(50, 1e5, 81.6), MyVector(), MyVector(.75, .75, .75), WOOD),//Botm
        Sphere(1e5, MyVector(50, -1e5 + 81.6, 81.6), MyVector(), MyVector(.75, .75, .75), DIFF),//Top
        Sphere(16.5, MyVector(27, 30, 47), MyVector(), MyVector(1, 1, 1) * .999, SPEC),//Mirr
        Sphere(16.5, MyVector(73, 35, 78), MyVector(), MyVector(1, 1, 1) * .999, REFR),//Glas
        Sphere(600, MyVector(50, 681.6 - .27, 81.6), MyVector(12, 12, 12), MyVector(), DIFF) //Lite
    };

#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP

    for (int y = 0; y < h; y++)
    { // Loop over image rows
        std::uniform_real_distribution<double> dis(0, 1);
        std::default_random_engine gen;
        auto ran = [&]
            {
                return dis(gen);
            };
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (int x = 0; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++ , r = MyVector())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * dis(gen), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * dis(gen), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        MyVector d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                            cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.dir;
                        r = r + radiance(Ray(cam.RayPos + d * 140, d.normalize()), 0, spheres, ran) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + MyVector(truncate(r.x), truncate(r.y), truncate(r.z)) * .25;
                }
    }
    std::cout << std::endl << clock() - begin << std::endl;


    Mat img(h, w, CV_8UC3, cv::Scalar(0, 0, 0));
    int rowNumber = img.rows;
    int colNumber = img.cols * img.channels();

    for (int i = 0; i < rowNumber; ++i)
    {
        uchar* data = img.ptr<uchar>(i);
        uchar* al = 0;
        uchar* ar = 0;
        for (int j = 0; j < colNumber; j = j + 3)
        {
            data[j] = toInt(c[i * w + j / 3].x);
            data[j + 1] = toInt(c[i * w + j / 3].y);
            data[j + 2] = toInt(c[i * w + j / 3].z);
        }
    }

    imwrite(".\\answer6.png", img);

    system("pause");
    return 0;
}

