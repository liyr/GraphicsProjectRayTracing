#pragma once
#include <cmath>
#include <ostream>

class MyVector
{
public:
    double x, y, z;

    explicit MyVector(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_)
    {
    }

    MyVector operator+(const MyVector& b) const
    {
        return MyVector(x + b.x, y + b.y, z + b.z);
    }

    MyVector operator-(const MyVector& b) const
    {
        return MyVector(x - b.x, y - b.y, z - b.z);
    }

    MyVector operator*(const double b) const
    {
        return MyVector(x * b, y * b, z * b);
    }

    MyVector operator*(const MyVector& b) const
    {
        return MyVector(x * b.x, y * b.y, z * b.z);
    }

    double norm() const
    {
        return (1 / sqrt(x * x + y * y + z * z));
    }

    double distance(const MyVector& b) const
    {
        return ((*this) - b).norm();
    }

    MyVector& normalize()
    {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }

    double operator | (const MyVector& b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }

    MyVector operator%(const MyVector& b) const
    {
        return MyVector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }

    void calSphereRepre(const MyVector& normal, const MyVector& bi_normal, double& theta, double& phi) const
    {
        MyVector tmp((*this)|(normal), (*this) | (bi_normal), (*this) | (normal % bi_normal));
        theta = acos(tmp.x);
        phi = atan2(tmp.z, tmp.y);
    }

    friend MyVector calSphereVec(const MyVector& normal, const MyVector& bi_normal, const double& theta, const double& phi)
    {
        return (normal * cos(theta) + bi_normal * sin(theta) * cos(phi) + (bi_normal % normal) * sin(theta) * sin(phi));
    }

    friend std::ostream& operator<<(std::ostream& out, MyVector& vec)
    {
        out << vec.x << " " << vec.y << " " << vec.z << std::endl;
        return out;
    }

    MyVector rotate(const MyVector& axis, const double angle) const
    {
        double temp;
        MyVector cross;
        double cos_ang = cos(angle);
        double sin_ang = sin(angle);


        MyVector out = (*this) * cos_ang;

        temp = operator|(axis);
        temp = temp * (1.0 - cos_ang);
        out = out + axis * temp;

        cross = axis % (*this);
        out = out + cross * sin_ang;
        return out;
    }

    bool inbox() const
    {
        auto tmp = 1e-4;
        return (0 - tmp < x && x < 100 + tmp && 0 - tmp < y && y < 81.6 + tmp && 0 - tmp < z && z < 170 + tmp);
    }
};

