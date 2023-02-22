#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "Matrix4.h"
#include "Vec3.h"

#include <string>

using namespace std;

class Camera
{

public:
    int cameraId;
    int projectionType; // 1 for perspective, 0 for orthographic
    Vec3 pos;
    Vec3 gaze;
    Vec3 u;
    Vec3 v;
    Vec3 w;
    double left, right, bottom, top;
    double near;
    double far;
    int horRes;
    int verRes;
    string outputFileName;

    Camera();

    Camera(int cameraId,
           int projectionType,
           Vec3 pos, Vec3 gaze,
           Vec3 u, Vec3 v, Vec3 w,
           double left, double right, double bottom, double top,
           double near, double far,
           int horRes, int verRes,
           string outputFileName);

    Camera(const Camera &other);

    friend std::ostream &operator<<(std::ostream &os, const Camera &c);

    Matrix4 get_camera_transformation_matrix();
    Matrix4 get_projection_transformation_matrix();
    Matrix4 get_viewport_transformation_matrix();
};

#endif
