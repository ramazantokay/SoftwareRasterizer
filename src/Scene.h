#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

	inline	double f_01(int x, int y, int x0, int x1, int y0, int y1)
	{
		return x * (y0 - y1) + y * (x1 - x0) + x0 * y1 - y0 * x1;
	}

	inline double f_12(int x, int y, int x1, int x2, int y1, int y2)
	{
		return x * (y1 - y2) + y * (x2 - x1) + x1 * y2 - y1 * x2;
	}

	inline double f_20(int x, int y, int x2, int x0, int y2, int y0)
	{
		return x * (y2 - y0) + y * (x0 - x2) + x2 * y0 - y2 * x0;
	}

public:
    Matrix4 geometric_transformation_matrix(int tranform_id, char type); // Rotating, Scaling, and Translating
	void line_rasterization(Vec4* v1, Vec4* v2, Camera& cam, Scene& scene);
	void triangle_rasterization(Vec4* v1, Vec4* v2, Vec4* v3, vector<Color> color, Camera& cam, Scene& scene);
	void liang_barsky_clipping(Vec4* v1, Vec4* v2);
	bool visible(double den, double num, double& t_e, double& t_l);
	void perspective_division(Vec4* v1, Vec4* v2, Vec4* v3);
	bool back_face_culling(Vec4* v1, Vec4* v2, Vec4* v3);
	Vec3 get_normal(Vec3 v1, Vec3 v2, Vec3 v3);

};

#endif
