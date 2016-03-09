#pragma once

#include "MyMath.h"

#ifndef MY_RAY
#define MY_RAY

struct Ray
{
	double t;
	Point3D p;	//position
	Point3D v;	//direct
};

#endif

#ifndef MY_OBJ
#define MY_OBJ
class Obj
{
public:
	virtual ~Obj(){}
};

class Triangle : public Obj
{
public:
	Point3D p0, p1, p2;	//point
	Point3D n;	//normal
	Point3D color;
	double d;
	double K_s;
	double K_t;

	void CalculateNormal(void);	//n = p1p0 x p1p2
	void SetPoints(Point3D p0, Point3D p1, Point3D p2);
	void SetColor(double r, double g, double b);
	void Calculate_d(void);	
	double f_of_Q(Point3D q);
	bool ContainmentTest(Point3D p);	//是否在三角型內, 是回傳1, 否回傳0

};

class Polygon : public Obj
{
public:
	Triangle t[2];

	void SetPoints(const Triangle &t0, const Triangle &t1);


};

#endif

Ray PrimeRay(const Point3D e, const Point3D p);
Ray ShadowRay(const Point3D p, const Point3D l);
Ray ReflectionRay(const Ray ray, const Point3D n, const Point3D p);
Ray TransmitRay(const Ray ray, const Point3D n, const Point3D p, double kt);