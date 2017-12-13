/*
assignment 2: implement ray tracing algorithm
name: Chenxi Wang
student id: 515030910467

the main function is at the bottom of this file
*/

#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>

using namespace std;

#define PI 3.141592653589793
#define INFINITY 1e8
#define MAX_RAY_DEPTH 5
float Image[640 * 480 * 4]; //used to display the data of the picture

// 3-dim vector template class
template<typename T>
class Vector3 {
public:
	T x, y, z; // coordinates
	Vector3() : x(T(0)), y(T(0)), z(T(0)) {}
	Vector3(T xx) : x(xx), y(xx), z(xx) {}
	Vector3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

	//vector normalization
	Vector3& normalize() {
		T nor2 = length2();
		if (nor2 > 0) {
			T invNor = 1 / sqrt(nor2);
			x *= invNor, y *= invNor, z *= invNor;
		}
		return *this;
	}

	//operator overload
	Vector3<T> operator * (const T &f) const { return Vector3<T>(x * f, y * f, z * f); }
	Vector3<T> operator * (const Vector3<T> &v) const { return Vector3<T>(x * v.x, y * v.y, z * v.z); }
	T dot(const Vector3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
	Vector3<T> outdot(const Vector3<T> &v) const { return Vector3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
	Vector3<T> operator - (const Vector3<T> &v) const { return Vector3<T>(x - v.x, y - v.y, z - v.z); }
	Vector3<T> operator + (const Vector3<T> &v) const { return Vector3<T>(x + v.x, y + v.y, z + v.z); }
	Vector3<T>& operator += (const Vector3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
	Vector3<T>& operator *= (const Vector3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	Vector3<T> operator - () const { return Vector3<T>(-x, -y, -z); }
	T length2() const { return x * x + y * y + z * z; }
	T length() const { return sqrt(length2()); }
	friend std::ostream & operator << (std::ostream &os, const Vector3<T> &v) {
		os << "[" << v.x << " " << v.y << " " << v.z << "]";
		return os;
	}
};
// class Vec3f is used to represent ver coordinate and dir vector
typedef Vector3<float> Vec3f;
// class Vec3d is used to represent the index of 3 vers in one mesh
typedef Vector3<int> Vec3d;

// an in-program-generated class Sphere
// a Sphere object can be reflective or diffuse
class Sphere {
public:
	Vec3f ct;								// the center coordinate of the sphere
	float r, r2;							// sphere radius and radius^2
	Vec3f surfColor, emisColor;				// surface color and emission light
	float transparency, reflection;         // surface transparency and reflectivity

	Sphere(const Vec3f &c, const float &rr, const Vec3f &sc,
		   const float &refl = 0, const float &transp = 0, const Vec3f &ec = 0):
		   ct(c), r(rr), r2(rr * rr), surfColor(sc), emisColor(ec),
		   transparency(transp), reflection(refl) {}

	// compute whether ray light has a intersection with the sphere
	// use 3d geometry rules
	// rs -> the position of the current ray source
	// rd -> the direction vector of the current ray
	// t0 and t1 are the two distances between the ray source and two intersections
	bool intersect(const Vec3f &rs, const Vec3f &rd, float &t0, float &t1) const {
		Vec3f l = ct - rs;
		float tca = l.dot(rd);
		// if the angle between the two vectors is greater than PI/2, there is no intersection
		if (tca < 0) return false;
		float d2 = l.dot(l) - tca * tca;
		// if the distance between sphere center and ray light is greater than r, no interseection
		if (d2 > r2) return false;
		float thc = sqrt(r2 - d2);
		// compute two distances and store them in t0 and t1
		t0 = tca - thc;
		t1 = tca + thc;

		return true;
	}
};


// a load-required class OFFModel
// load model from .off file(Priceton object format)
// no reflectivity in OFFModel
// an OFFModel object has several vers, faces and corresponding face norms
// the method to judge intersection is different from Sphere
class OFFModel {
public:
	int nVer, nFace;				// the number of vers and faces
	vector<Vec3f> vers;				// a series of vers contained in model
	vector<Vec3d> faces;			// a series of ver idxs to construct faces, each face containing 3 vers
	vector<Vec3f> faceNorms;		// the norms of the corresponding faces
	Vec3f surfColor;				// the surface color of this model

	// use .off file to generate model
	// the format of .off file is easy to parse
	OFFModel(char* filename) {
		ifstream fin(filename);
		char tmpchar[100];
		int tmpnum;
		// define the model surface color
		surfColor.x = 1;  surfColor.y = 0; surfColor.z = 0;
		fin >> tmpchar;
		fin >> nVer >> nFace >> tmpnum;
		// get the ver information
		for (int i = 0; i < nVer; ++i) {
			Vec3f pt;
			fin >> pt.x >> pt.y >> pt.z;
			pt.z = -pt.z;
			pt = pt*0.13 + Vec3f(-3, -6, -40);
			vers.push_back(pt);
		}
		// get the face information (including face norms)
		for (int i = 0; i < nFace; ++i) {
			Vec3d fs;
			Vec3f nm;
			fin >> tmpnum >> fs.x >> fs.y >> fs.z;
			faces.push_back(fs);
			// compute the face norm using outter production
			Vec3f v1 = vers[fs.y] - vers[fs.x], v2 = vers[fs.z] - vers[fs.y];
			nm = v1.outdot(v2);
			nm.normalize();
			faceNorms.push_back(nm);
		}
		fin.close();
	}

	// public function for OFFModel to judge and compute the intersection
	// each time running this function, it will travers all faces to get the result
	int intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t) {
		// judge and compute the intersection for one single face
		for (int i = 0; i < nFace; ++i) {
			if (intersect(i, rayorig, raydir, t))
				return i;
		}
		// if no intersection, return a negetive number -1
		return -1;
	}

private:
	// private function for a single face to judge and compute the intersection
	// first get the intersection in the face plane, then judge if it's in the face triangle
	// i -> the current face index
	// rs -> the position of the current ray source
	// rd -> the direction vector of the current ray
	bool intersect(int i, const Vec3f &rs, const Vec3f &rd, float &t) {
		// compute the intersection of the current face plane and ray light
		Vec3f nm = faceNorms[i], pt1 = vers[faces[i].x], pt2 = vers[faces[i].y], pt3 = vers[faces[i].z];
		float tmp_t;
		float dir_dot_norm = rd.dot(nm);
		if (dir_dot_norm >= 0) {
			return false;
		}
		tmp_t = nm.dot(pt1 - rs) * (1 / dir_dot_norm);

		Vec3f P = rs + rd * tmp_t;
		// if P is in face, the traingle area should be
		// equal to the sum of 3 sub traingles containing P and two vers in the face
		float areaAll = ((pt2 - pt1).outdot(pt3 - pt2)).length();
		float area1 = ((pt2 - pt1).outdot(P - pt2)).length();
		float area2 = ((pt3 - pt2).outdot(P - pt3)).length();
		float area3 = ((pt1 - pt3).outdot(P - pt1)).length();
		if (area1 + area2 + area3 > areaAll + area1 / 1000) {
			// P is not in the triangle
			return false;
		}

		t = tmp_t;
		return true;
	}
};

// mix the colors
float mix(const float &a, const float &b, const float &mix) {
	return b * mix + a * (1 - mix);
}


// ray tracing implementation
// different processing method for sphere and model
// a recursive function
Vec3f trace(const Vec3f &rs, const Vec3f &rd, const std::vector<Sphere> &spheres, OFFModel &model, const int &depth) {
	float tnear = INFINITY;
	const Sphere* sphere = NULL;
	// find intersection of this ray and spheres
	for (unsigned i = 0; i < spheres.size(); ++i) {
		float t0 = INFINITY, t1 = INFINITY;
		if (spheres[i].intersect(rs, rd, t0, t1)) {
			if (t0 < 0) t0 = t1;
			if (t0 < tnear) {
				tnear = t0;
				sphere = &spheres[i];
			}
		}
	}
	// compute the intersection of the current ray and model surface
	float t = INFINITY;
	int face_id = model.intersect(rs, rd, t);
	if (t < tnear) {
		tnear = t;
		sphere = NULL;
	}

	Vec3f surfColor = 0;

	// if intersction belongs to model, return the light from model
	if (face_id >= 0) {
		surfColor = model.surfColor * model.faceNorms[face_id].dot(-rd);
		return surfColor;
	}
	// if no intersection, return background color
	if (!sphere) { return Vec3f(2); }

	// if intersection belongs to sphere, recursively compute the color
	Vec3f P = rs + rd * tnear;
	Vec3f N = P - sphere->ct;
	N.normalize();
	float bias = 1e-4; // just a bias for comoputing
	bool inside = false;
	if (rs.dot(N) > 0) { N = -N; inside = true; }
	if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
		float facingratio = -rd.dot(N);
		float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
		// compute reflection direction
		Vec3f refldir = rd - N * 2 * rd.dot(N);
		refldir.normalize();
		Vec3f reflection = trace(P + N * bias, refldir, spheres, model, depth + 1);
		Vec3f refraction = 0;
		// if the sphere is also transparent, compute refraction ray
		if (sphere->transparency) {
			float ior = 1.1, eta = (inside) ? ior : 1 / ior;
			float cosi = -N.dot(rd);
			float k = 1 - eta * eta * (1 - cosi * cosi);
			Vec3f refrdir = rd * eta + N * (eta *  cosi - sqrt(k));
			refrdir.normalize();
			refraction = trace(P - N * bias, refrdir, spheres, model, depth + 1);
		}
		// the result is a mix of reflection and refraction (if the sphere is transparent)
		surfColor = (reflection*fresneleffect + refraction*(1 - fresneleffect)*sphere->transparency)*sphere->surfColor;
	}
	else {
		// the sphere is diffuse, we don't have any more recursive computation
		for (unsigned i = 0; i < spheres.size(); ++i) {
			if (spheres[i].emisColor.x > 0) {
				// this is a light, we need to judge if there is any shading
				Vec3f transmission = 1;
				Vec3f lightDirection = spheres[i].ct - P;
				lightDirection.normalize();
				// from light source to judge the intersection with sphere
				for (unsigned j = 0; j < spheres.size(); ++j) {
					if (i != j) {
						float t0, t1;
						if (spheres[j].intersect(P + N * bias, lightDirection, t0, t1)) {
							transmission = 0;
							break;
						}
					}
				}
				// from light source to judge the intersection with model
				float tmp_t;
				int tmp_face_id = model.intersect(P + N * bias, lightDirection, tmp_t);
				if (tmp_face_id > 0) transmission = 0;
				surfColor += sphere->surfColor * transmission * max(float(0), N.dot(lightDirection)) * spheres[i].emisColor;
			}
		}
	}

	return surfColor + sphere->emisColor;
}


// render the scene and get the color of each pixel
void render(const std::vector<Sphere> &spheres, OFFModel model) {
	unsigned width = 640, height = 480;
	Vec3f *image = new Vec3f[width * height], *pixel = image;
	float invWidth = 1 / float(width), invHeight = 1 / float(height);
	float fov = 30, aspectratio = width / float(height);
	float angle = tan(PI * 0.5 * fov / 180.);

	for (unsigned y = 0; y < height; ++y) {
		for (unsigned x = 0; x < width; ++x, ++pixel) {
			float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
			float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
			Vec3f raydir(xx, yy, -1);
			raydir.normalize();
			// trace rays
			*pixel = trace(Vec3f(0, 0, 0), raydir, spheres, model, 0);
		}
	}

	//store the data in the Image array
	float *LinePD = Image;
	for (int i = height - 1; i >= 0; --i) {
		for (int j = 0; j < width; ++j) {
			LinePD[0] = min(float(1), image[j + i * width].x);
			LinePD[1] = min(float(1), image[j + i * width].y);
			LinePD[2] = min(float(1), image[j + i * width].z);
			LinePD[3] = 0.02;
			LinePD += 4;
		}
	}
	delete[] image;
}

//display the rendered image using OpenGL
void display() {
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(640, 480, GL_RGBA, GL_FLOAT, Image);
	glFlush();
}


int main(int argc, char **argv) {

	std::vector<Sphere> spheres;
	// construct spheres
	// use a huge diffuse sphere to represent the plane
	spheres.push_back(Sphere(Vec3f(0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
	// construct two reflective spheres
	spheres.push_back(Sphere(Vec3f(7.0, -1, -20), 2, Vec3f(1.0, 0.76, 0.0), 1, 0.0));
	spheres.push_back(Sphere(Vec3f(8.0, 0, -30), 3, Vec3f(0., 1, 0.), 1, 0.0));
	// load .off model from the exact file
	OFFModel bunny("model/armadillo_32.off");
	// used as a light
	spheres.push_back(Sphere(Vec3f(0.0, 20, -30), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
	// render the scene
	render(spheres, bunny);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowSize(640, 480);
	glutInitWindowPosition(200, 200);
	glutCreateWindow("HW2-ChenxiWang");
	glClearColor(1, 1, 1, 1);
	glutDisplayFunc(display);
	glutMainLoop();

	return 0;
}
