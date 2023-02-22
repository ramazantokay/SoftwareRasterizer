#include "Camera.h"
#include "Helpers.h"
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

Camera::Camera() {}

Camera::Camera(int cameraId,
               int projectionType,
               Vec3 pos, Vec3 gaze,
               Vec3 u, Vec3 v, Vec3 w,
               double left, double right, double bottom, double top,
               double near, double far,
               int horRes, int verRes,
               string outputFileName)
{

    this->cameraId = cameraId;
    this->projectionType = projectionType;
    this->pos = pos;
    this->gaze = gaze;
    this->u = u;
    this->v = v;
    this->w = w;
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->near = near;
    this->far = far;
    this->horRes = horRes;
    this->verRes = verRes;
    this->outputFileName = outputFileName;
}

Camera::Camera(const Camera &other)
{
    this->cameraId = other.cameraId;
    this->projectionType = other.projectionType;
    this->pos = other.pos;
    this->gaze = other.gaze;
    this->u = other.u;
    this->v = other.v;
    this->w = other.w;
    this->left = other.left;
    this->right = other.right;
    this->bottom = other.bottom;
    this->top = other.top;
    this->near = other.near;
    this->far = other.far;
    this->horRes = other.horRes;
    this->verRes = other.verRes;
    this->outputFileName = other.outputFileName;
}

ostream &operator<<(ostream &os, const Camera &c)
{
    const char *camType = c.projectionType ? "perspective" : "orthographic";

    os << fixed << setprecision(6) << "Camera " << c.cameraId << " (" << camType << ") => pos: " << c.pos << " gaze: " << c.gaze << endl
       << "\tu: " << c.u << " v: " << c.v << " w: " << c.w << endl
       << fixed << setprecision(3) << "\tleft: " << c.left << " right: " << c.right << " bottom: " << c.bottom << " top: " << c.top << endl
       << "\tnear: " << c.near << " far: " << c.far << " resolutions: " << c.horRes << "x" << c.verRes << " fileName: " << c.outputFileName;

    return os;
}

Matrix4 Camera::get_camera_transformation_matrix() // M_cam matrix (7.4) equation from the book
{
    double m_camera_transformation_matrix[4][4] = {
        {u.x, u.y, u.z, 0},
        {v.x, v.y, v.z, 0},
        {w.x, w.y, w.z, 0},
        {0, 0, 0, 1},
    };

    double m_aligned_e[4][4] = {
        {1, 0, 0, -pos.x},
        {0, 1, 0, -pos.y},
        {0, 0, 1, -pos.z},
        {0, 0, 0, 1},
    };

    Matrix4 m_1(m_camera_transformation_matrix);
    Matrix4 m_2(m_aligned_e);
    Matrix4 m_cam_tran_matrix = getIdentityMatrix();

    m_cam_tran_matrix = multiplyMatrixWithMatrix(m_2, m_cam_tran_matrix);
    m_cam_tran_matrix = multiplyMatrixWithMatrix(m_1, m_cam_tran_matrix);

    return m_cam_tran_matrix;
}

Matrix4 Camera::get_projection_transformation_matrix()
{
    double abs_near = abs(near);
    double abs_far = abs(far);
    double one_div_abs_near_far = 1 / (abs_near - abs_far);

    double one_div_top_bottom = 1 / (top - bottom);
    double one_div_far_near = 1 / (far - near);
    double one_div_right_left = 1 / (right - left);

    if (projectionType == 0) // for orthografic, from the book page 143
    {
        double d_orthographic[4][4] =
            {
                {2 * one_div_right_left, 0, 0, -((right + left) * one_div_right_left)},
                {0, 2 * one_div_top_bottom, 0, -((top + bottom) * one_div_top_bottom)},
                {0, 0, -2 * one_div_far_near, ((near + far) * one_div_far_near)},
                {0, 0, 0, 1}};
        Matrix4 m_orthographic(d_orthographic);
        return m_orthographic;
    }
    else // for perspective, from the book page 153
    {
        double d_perspective[4][4] =
            {
                {2 * abs_near * one_div_right_left, 0, (right + left) * one_div_right_left, 0},
                {0, 2 * abs_near * one_div_top_bottom, (top + bottom) * one_div_top_bottom, 0},
                {0, 0, (abs_near + abs_far) * one_div_abs_near_far, 2 * abs_far * abs_near * one_div_abs_near_far},
                {0, 0, -1, 0}};

        Matrix4 m_perspective(d_perspective);
        return m_perspective;
    }
}

Matrix4 Camera::get_viewport_transformation_matrix()
{
    double d_viewport[4][4] = {
        {horRes * 0.5, 0, 0, (horRes - 1) * 0.5},
        {0, verRes * 0.5, 0, (verRes - 1) * 0.5},
        {0, 0, 1, 0},
        {0, 0, 0, 1}};
    Matrix4 m_viewport(d_viewport);
    return m_viewport;
}
