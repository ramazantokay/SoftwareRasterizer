#ifndef __COLOR_H__
#define __COLOR_H__
#include "Vec3.h"

#include <iostream>
class Vec3;

class Color
{
public:
    double r, g, b;

    Color();
    Color(double r, double g, double b);
    Color(const Color &other);
    Color(const Vec3 &other);
    friend std::ostream& operator<<(std::ostream& os, const Color& c);
};

#endif
