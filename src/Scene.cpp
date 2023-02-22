#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/

Matrix4 Scene::geometric_transformation_matrix(int tranform_id, char type)
{
	// cout << "yooovvv geo transofmasyonnn" << endl;
	// cout << "¯\\_(ツ)_/¯" << endl;
	if (type == 't') // Translation
	{
		Translation *m_translation = translations[tranform_id];
		double d_t[4][4] = {
			{1, 0, 0, m_translation->tx},
			{0, 1, 0, m_translation->ty},
			{0, 0, 1, m_translation->tz},
			{0, 0, 0, 1}};
		Matrix4 m_t(d_t);
		return m_t;
	}
	else if (type == 's') // Scaling
	{
		Scaling *m_scaling = scalings[tranform_id];
		double d_s[4][4] =
			{
				{m_scaling->sx, 0, 0, 0},
				{0, m_scaling->sy, 0, 0},
				{0, 0, m_scaling->sz, 0},
				{0, 0, 0, 1},
			};
		Matrix4 s_t(d_s);
		return s_t;
	}
	else if (type == 'r') // Rotating
	{
		Rotation *m_rotation = rotations[tranform_id];
		double angle = (m_rotation->angle / 180) * M_PI;
		double cos_angle = cos(angle);
		double sin_angle = sin(angle);

		//@TODO: continue to implment rotation.

		// Find v, set the smallest element of u in a absolute sense
		Vec3 v;
		// Alternative method from Modelling Transformaiton slide page 41.
		Vec3 u = Vec3(m_rotation->ux, m_rotation->uy, m_rotation->uz, 0);
		Vec3 abs_u = Vec3(abs(m_rotation->ux), abs(m_rotation->uy), abs(m_rotation->uz), 0);

		// Find smallest element of u
		if (abs_u.x <= abs_u.y && abs_u.x <= abs_u.z) // For the x
		{
			v.x = 0;
			v.y = -u.z;
			v.z = u.y;
		}
		else if (abs_u.y <= abs_u.x && abs_u.y <= abs_u.z) // For the y
		{
			v.y = 0;
			v.x = -u.z;
			v.z = u.x;
		}
		else if (abs_u.z <= abs_u.x && abs_u.z <= abs_u.y) // For the z
		{
			v.z = 0;
			v.x = -u.y;
			v.y = u.x;
		}

		// Find w by cross product of u and v
		Vec3 w = crossProductVec3(u, v);
		v = normalizeVec3(v);
		w = normalizeVec3(w);

		// inverse orthonormal matrix
		double d_inverse_orth_matrix[4][4] = {
			{u.x, v.x, w.x, 0},
			{u.y, v.y, w.y, 0},
			{u.z, v.z, w.z, 0},
			{0, 0, 0, 1},
		};
		// rotation x matrix
		double d_rotation_x_matrix[4][4] = {
			{1, 0, 0, 0},
			{0, cos_angle, -sin_angle, 0},
			{0, sin_angle, cos_angle, 0},
			{0, 0, 0, 1},
		};
		// orthonormal matrix
		double d_orth_matrix[4][4] = {
			{u.x, u.y, u.z, 0},
			{v.x, v.y, v.z, 0},
			{w.x, w.y, w.z, 0},
			{0, 0, 0, 1},
		};

		Matrix4 m_inverse_orth(d_inverse_orth_matrix);
		Matrix4 m_orth(d_orth_matrix);
		Matrix4 m_rotation_x(d_rotation_x_matrix);

		Matrix4 m_final_rotation = getIdentityMatrix();
		m_final_rotation = multiplyMatrixWithMatrix(m_orth, m_final_rotation);
		m_final_rotation = multiplyMatrixWithMatrix(m_rotation_x, m_final_rotation);
		m_final_rotation = multiplyMatrixWithMatrix(m_inverse_orth, m_final_rotation);

		// @TODO create debug printers
		// cout << "¯(ㆆ _ ㆆ)" << endl;

		// cout << "(̿▀̿ ̿Ĺ̯̿̿▀̿ ̿)̄" << endl;

		return m_final_rotation;
	}
	else
		throw std::runtime_error("baba bi sorun vaar ama bi kontrol et (̿▀̿ ̿Ĺ̯̿̿▀̿ ̿)̄");
}

void Scene::triangle_rasterization(Vec4 *v1, Vec4 *v2, Vec4 *v3, vector<Color> color, Camera &cam, Scene &scene)
{
	double x_min, x_max, y_min, y_max;

	x_min = min(round(v1->x), min(round(v2->x), round(v3->x)));
	x_max = max(round(v1->x), max(round(v2->x), round(v3->x)));

	y_min = min(round(v1->y), min(round(v2->y), round(v3->y)));
	y_max = max(round(v1->y), max(round(v2->y), round(v3->y)));

	for (int y = y_min; y < y_max; y++)
	{
		for (int x = x_min; x < x_max; x++)
		{
			if (x < 0 || x >= cam.horRes || y < 0 || y >= cam.verRes)
				continue;

			double alpha, beta, gamma;

			alpha = f_12(x, y, v2->x, v3->x, v2->y, v3->y) / f_12(v1->x, v1->y, v2->x, v3->x, v2->y, v3->y);
			beta = f_20(x, y, v3->x, v1->x, v3->y, v1->y) / f_20(v2->x, v2->y, v3->x, v1->x, v3->y, v1->y);
			gamma = f_01(x, y, v1->x, v2->x, v1->y, v2->y) / f_01(v3->x, v3->y, v1->x, v2->x, v1->y, v2->y);

			if (alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				Color temp_color;
				temp_color.r = round(alpha * color[0].r + beta * color[1].r + gamma * color[2].r);
				temp_color.g = round(alpha * color[0].g + beta * color[1].g + gamma * color[2].g);
				temp_color.b = round(alpha * color[0].b + beta * color[1].b + gamma * color[2].b);
				scene.image[x][y] = temp_color;
			}
		}
	}
}

void Scene::line_rasterization(Vec4 *v1, Vec4 *v2, Camera &cam, Scene &scene)
{
	Vec4 *t_1 = (v1->x <= v2->x) ? v1 : v2;
	Vec4 *t_2 = (v1->x > v2->x) ? v1 : v2;

	Color color_0 = *scene.colorsOfVertices[t_1->colorId - 1];
	Color color_1 = *scene.colorsOfVertices[t_2->colorId - 1];

	Vec3 c0(color_0);
	Vec3 c1(color_1);
	Vec3 c;

	// cout << "nerde seg fault yiyoruma acaba line_raster" << endl;

	double x_0 = t_1->x;
	double x_1 = t_2->x;

	double y_0 = t_1->y;
	double y_1 = t_2->y;

	double slope = (y_1 - y_0) / (x_1 - x_0);

	// cout << "nerde seg fault yiyoruma acaba line_raster slope" << endl;

	if (slope >= 0 && slope <= 1)
	{
		double y = y_0;
		double d = 2 * (y_0 - y_1) + (x_1 - x_0);
		c = c0;
		Vec3 cd = multiplyVec3WithScalar(subtractVec3(c1, c0), 1 / (x_1 - x_0));

		for (int x = (int)x_0; x <= (int)x_1; x++)
		{

			Color temp_color(c);
			// cout << "nerde seg fault yiyoruma acaba line_raster slope >=0 1" << endl;

			if (x >= 0 && x < image.size() && y >= 0 && y < scene.image[0].size()) // bu ekleyince seg fault gitti, vay arkadas :|
				scene.image[x][y] = temp_color;
			// cout << "nerde seg fault yiyoruma acaba line_raster slope >=0 2" << endl;

			if (d < 0)
			{
				y++;
				d += 2 * ((y_0 - y_1) + (x_1 - x_0));
			}
			else
			{
				d += 2 * ((y_0 - y_1));
			}
			c = addVec3(c, cd);
		}
	}
	else if (slope > 1)
	{
		double x = x_0;
		double d = 2 * (x_0 - x_1) + (y_1 - y_0);

		c = c0;
		Vec3 cd = multiplyVec3WithScalar(subtractVec3(c1, c0), 1 / (y_1 - y_0));

		for (int y = (int)y_0; y <= (int)y_1; y++)
		{
			Color temp_color(c);
			// cout << "nerde seg fault yiyoruma acaba line_raster slope>1 222" << endl;
			if (x >= 0 && x < image.size() && y >= 0 && y < scene.image[0].size())
				scene.image[x][y] = temp_color;

			// cout << "nerde seg fault yiyoruma acaba line_raster slope>1" << endl;

			if (d < 0)
			{
				x++;
				d += 2 * ((x_0 - x_1) + (y_1 - y_0));
			}
			else
			{
				d += 2 * ((x_0 - x_1));
			}
			c = addVec3(c, cd);
		}
	}
	else if (slope < 0 && slope >= -1)
	{
		double y = y_1;
		double d = 2 * (y_1 - y_0) + (x_1 - x_0);

		c = c1;
		Vec3 cd = multiplyVec3WithScalar(subtractVec3(c0, c1), 1 / (x_1 - x_0));

		for (int x = (int)x_1; x >= (int)x_0; x--)
		{
			Color temp_color(c);

			// cout << "nerde seg fault yiyoruma acaba line_raster slope<0 >=-1 1" << endl;
			if (x >= 0 && x < image.size() && y >= 0 && y < scene.image[0].size())
				scene.image[x][y] = temp_color;
			// cout << "nerde seg fault yiyoruma acaba line_raster slope<0 >=-1 2" << endl;

			if (d < 0)
			{
				y++;
				d += 2 * ((y_1 - y_0 + x_1 - x_0));
			}
			else
			{
				d += 2 * (y_1 - y_0);
			}
			c = addVec3(c, cd);
		}
	}
	else if (slope < -1)
	{
		double x = x_1;
		double d = 2 * (x_0 - x_1) + (y_0 - y_1);

		c = c1;
		Vec3 cd = multiplyVec3WithScalar(subtractVec3(c0, c1), 1 / (y_0 - y_1));

		for (int y = (int)y_1; y <= (int)y_0; y++)
		{
			Color temp_color(c);
			// cout << "nerde seg fault yiyoruma acaba line_raster slope<-1 1" << endl;

			if (x >= 0 && x < image.size() && y >= 0 && y < scene.image[0].size())
				scene.image[x][y] = temp_color;

			// cout << "nerde seg fault yiyoruma acaba line_raster slope<-1 2" << endl;

			if (d < 0)
			{
				x--;
				d += 2 * (x_0 - x_1 + y_0 - y_1);
			}
			else
			{
				d += 2 * ((x_0 - x_1));
			}
			c = addVec3(c, cd);
		}
	}
}


bool Scene::visible(double den, double num, double &t_e, double &t_l)
{
	double d_t;
	if (den > 0) // potentially entering
	{
		d_t = num / den;
		if (d_t > t_l)
			return false;
		if (d_t > t_e)
			t_e = d_t;
	}
	else if (den < 0) // potentially leaving
	{
		d_t = num / den;
		if (d_t < t_e)
			return false;
		if (d_t < t_l)
			t_l = d_t;
	}
	else if (num > 0) // line paralell to edge
	{
		return false;
	}
	return true;
}

void Scene::liang_barsky_clipping(Vec4 *v1, Vec4 *v2)
{
	double t_e = 0;
	double t_l = 1;
	Vec3 d((v2->x - v1->x), (v2->y - v1->y), (v2->z - v1->z), 0);

	if (visible(d.x, -1 - v1->x, t_e, t_l))
	{
		cout << "1.if "<< endl;
		if (visible(-d.x, v1->x - 1, t_e, t_l))
		{
		cout << "2.if "<< endl;

			if (visible(d.y, -1 - v1->y, t_e, t_l))
			{
			cout << "3.if "<< endl;

				if (visible(-d.y, v1->y - 1, t_e, t_l))
				{
				cout << "4.if "<< endl;

					if (visible(d.z, -1 - v1->z, t_e, t_l))
					{
					cout << "5.if "<< endl;

						if (visible(-d.z, v1->z - 1, t_e, t_l))
						{
							cout << "6.if "<< endl;
							if (t_l < 1)
							{
								cout << "t_l.if "<< endl;

								v2->x = v1->x + d.x * t_l;
								v2->y = v1->y + d.y * t_l;
								v2->z = v1->z + d.z * t_l;
							}
							if (t_e > 0)
							{
								cout << "t_e.if "<< endl;

								v1->x = v2->x + d.x * t_e;
								v1->y = v2->y + d.y * t_e;
								v1->z = v2->z + d.z * t_e;
							}
						}
					}
				}
			}
		}
	}
}

void Scene::perspective_division(Vec4 *v1, Vec4 *v2, Vec4 *v3)
{
	v1->x = v1->x / v1->t;
	v1->y = v1->y / v1->t;
	v1->z = v1->z / v1->t;
	v1->t = 1;

	v2->x = v2->x / v2->t;
	v2->y = v2->y / v2->t;
	v2->z = v2->z / v2->t;
	v2->t = 1;

	v3->x = v3->x / v3->t;
	v3->y = v3->y / v3->t;
	v3->z = v3->z / v3->t;
	v3->t = 1;
}

Vec3 Scene::get_normal(Vec3 v1, Vec3 v2, Vec3 v3)
{
	Vec3 u = subtractVec3(v2, v1);
	Vec3 v = subtractVec3(v3, v1);
	Vec3 n = crossProductVec3(u, v);
	n = normalizeVec3(n);
	return n;
}

bool Scene::back_face_culling(Vec4 *v1, Vec4 *v2, Vec4 *v3)
{
	Vec3 bfc_1(v1->x, v1->y, v1->z, 0);
	Vec3 bfc_2(v2->x, v2->y, v2->z, 0);
	Vec3 bfc_3(v3->x, v3->y, v3->z, 0);

	Vec3 n = get_normal(bfc_1, bfc_2, bfc_3);
	Vec3 view = subtractVec3(bfc_1, Vec3(0, 0, 0, 0)); // from 7th slide page 66

	return (dotProductVec3(n, view) < 0) ? true : false;
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	for (Mesh *mesh : meshes)
	{
		for (Triangle &triangle : mesh->triangles)
		{
			Vec3 *vertex_1, *vertex_2, *vertex_3;
			vertex_1 = vertices[triangle.vertexIds[0] - 1];
			vertex_2 = vertices[triangle.vertexIds[1] - 1];
			vertex_3 = vertices[triangle.vertexIds[2] - 1];

			// cout << "Vertcies 3" << endl;
			// cout << "v1: " << *vertex_1 << " v2: " << *vertex_2 << " v3: " << *vertex_3 << endl;

			Vec4 *ver_1, *ver_2, *ver_3;
			ver_1 = new Vec4(vertex_1->x, vertex_1->y, vertex_1->z, 1, vertex_1->colorId);
			ver_2 = new Vec4(vertex_2->x, vertex_2->y, vertex_2->z, 1, vertex_2->colorId);
			ver_3 = new Vec4(vertex_3->x, vertex_3->y, vertex_3->z, 1, vertex_3->colorId);

			// cout << "Before Modelling Transformation" << endl;
			// cout << "v1: " << *ver_1 << " v2: " << *ver_2 << " v3: " << *ver_3 << endl;

			// ******************** Modelling Transformation Begins ********************
			Matrix4 m_general_matrix = getIdentityMatrix();
			for (int i = 0; i < mesh->numberOfTransformations; i++)
			{
				m_general_matrix = multiplyMatrixWithMatrix(geometric_transformation_matrix(mesh->transformationIds[i] - 1, mesh->transformationTypes[i]), m_general_matrix);
			}

			*ver_1 = multiplyMatrixWithVec4(m_general_matrix, *ver_1);
			*ver_2 = multiplyMatrixWithVec4(m_general_matrix, *ver_2);
			*ver_3 = multiplyMatrixWithVec4(m_general_matrix, *ver_3);

			cout << "After Modelling Transformation" << endl;
			cout << "ver_1: " << *ver_1 << " ver_2: " << *ver_2 << " ver_3: " << *ver_3 << endl;
			cout << "¯(ㆆ _ ㆆ)" << endl;
			// ******************** Modelling Transformation End ********************

			// ******************** Camera Tranformation Begin ********************

			Matrix4 m_camera_trans_matrix = camera->get_camera_transformation_matrix();
			m_general_matrix = multiplyMatrixWithMatrix(m_camera_trans_matrix, m_general_matrix);

			*ver_1 = multiplyMatrixWithVec4(m_camera_trans_matrix, *ver_1);
			*ver_2 = multiplyMatrixWithVec4(m_camera_trans_matrix, *ver_2);
			*ver_3 = multiplyMatrixWithVec4(m_camera_trans_matrix, *ver_3);

			cout << "After Camera Transformation" << endl;
			cout << "ver_1: " << *ver_1 << " ver_2: " << *ver_2 << " ver_3: " << *ver_3 << endl;
			cout << "¯(ㆆ _ ㆆ)" << endl;

			// ******************** Camera Tranformation End ********************

			// ******************** Projection Tranformation Begin ********************

			Matrix4 m_projection_matrix = camera->get_projection_transformation_matrix();
			m_general_matrix = multiplyMatrixWithMatrix(m_projection_matrix, m_general_matrix);

			*ver_1 = multiplyMatrixWithVec4(m_projection_matrix, *ver_1);
			*ver_2 = multiplyMatrixWithVec4(m_projection_matrix, *ver_2);
			*ver_3 = multiplyMatrixWithVec4(m_projection_matrix, *ver_3);

			cout << "After Projection Transformation" << endl;
			cout << "ver_1: " << *ver_1 << " ver_2: " << *ver_2 << " ver_3: " << *ver_3 << endl;
			cout << "¯(ㆆ _ ㆆ)" << endl;

			// ******************** Projection Tranformation End ********************

			// ******************** Backface Culling Begin ********************

			if (cullingEnabled)
			{
				if (back_face_culling(ver_1, ver_2, ver_3))
				{
					cout << "Culling is done" << endl;
					continue;
				}
			}
			// ******************** Backface Culling End ********************

			// ******************** Rasterize Begin ********************

			if (mesh->type == 1) // Solid
			{
				// perspective division
				perspective_division(ver_1, ver_2, ver_3);

				// viewport transformation
				Matrix4 m_viewport_matrix = camera->get_viewport_transformation_matrix();
				m_general_matrix = multiplyMatrixWithMatrix(m_viewport_matrix, m_general_matrix);

				*ver_1 = multiplyMatrixWithVec4(m_viewport_matrix, *ver_1);
				*ver_2 = multiplyMatrixWithVec4(m_viewport_matrix, *ver_2);
				*ver_3 = multiplyMatrixWithVec4(m_viewport_matrix, *ver_3);

				// get color
				vector<Color> color_vector;
				color_vector.push_back(*colorsOfVertices[ver_1->colorId - 1]);
				color_vector.push_back(*colorsOfVertices[ver_2->colorId - 1]);
				color_vector.push_back(*colorsOfVertices[ver_3->colorId - 1]);

				triangle_rasterization(ver_1, ver_2, ver_3, color_vector, *camera, *this);
			}
			else // Wireframe
			{
				// @TODO: wireframe implementation


				// cout << "nerde seg fault yiyoruma acaba 1" << endl;
				// clipping
				// clipping algosu doğru ama neden saçmalıyor anlamış değilim valla :|
				liang_barsky_clipping(ver_1, ver_2);
				liang_barsky_clipping(ver_2, ver_3);
				liang_barsky_clipping(ver_3, ver_1);

				// perspective division
				perspective_division(ver_1, ver_2, ver_3);
				
				// viewport transformation
				Matrix4 m_viewport_matrix = camera->get_viewport_transformation_matrix();
				m_general_matrix = multiplyMatrixWithMatrix(m_viewport_matrix, m_general_matrix);

				*ver_1 = multiplyMatrixWithVec4(m_viewport_matrix, *ver_1);
				*ver_2 = multiplyMatrixWithVec4(m_viewport_matrix, *ver_2);
				*ver_3 = multiplyMatrixWithVec4(m_viewport_matrix, *ver_3);

				// cout << "nerde seg fault yiyoruma acaba 2" << endl;

				line_rasterization(ver_1, ver_2, *camera, *this);
				line_rasterization(ver_2, ver_3, *camera, *this);
				line_rasterization(ver_3, ver_1, *camera, *this);

				// cout << "nerde seg fault yiyoruma acaba 6" << endl;
			}

			// ******************** Rasterize End ********************
		}
	}
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL)
	{
		str = pElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			cullingEnabled = true;
		}
		else
		{
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			cam->projectionType = 0;
		}
		else
		{
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = 0;
		}
		else
		{
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
		str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}
