#pragma once

#ifndef FIRST_DEFINE
#define FIRST_DEFINE
#include <vector>
#include <cstdlib>
#include <iostream>
#include <GL/glut.h>
#include <cmath>
using namespace std;

#endif

#define PI 3.14159
#define ANGLE PI/180

#ifndef P_3D
#define P_3D
struct Point3D
{
    double x, y, z;

};
Point3D operator-(const Point3D &a, const Point3D &b);	//return a-b
Point3D operator-(const double v, const Point3D &b);	//return v-b
Point3D operator+(const Point3D &a, const Point3D &b);	//return a+b
Point3D operator*(const double v, const Point3D &b);	//return v*b
Point3D operator*(const Point3D &a, const Point3D &b);	//return a*b


#endif
/*
void Swap(double &a, double &b)
{
	double k; 
	k = a ; a = b; b = k;
}*/

void NormalizeV(Point3D &t);
double InnerProduct(const Point3D &v1, const Point3D &v2);
void CrossProduct(Point3D &u, const Point3D &w, const Point3D &view_up);	//w x view_up = u
double CalcualtionLength(const Point3D &p0, const Point3D &p1);	// | p0-p1 |
Point3D DiffuseColor(Point3D color, const Point3D l, const Point3D n);
Point3D SpecularColor(Point3D color, const Point3D v, const Point3D l, const Point3D n);


//void Gauss(double *value, double *cons, int dimension);
