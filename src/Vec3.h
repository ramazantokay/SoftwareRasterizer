#ifndef __VEC3_H__
#define __VEC3_H__

#include "Color.h" 
#include <iostream>
using namespace std;


class Color;

class Vec3
{
public:
    double x, y, z;
    int colorId;

    Vec3();
    Vec3(double x, double y, double z, int colorId);
    Vec3(const Vec3 &other);
    
    Vec3(const Color &other);

    double getElementAt(int index);
    
    friend std::ostream& operator<<(std::ostream& os, const Vec3& v);
};

#endif