#include "MyTool.h"

void DisplayFunc(void);
void RayTrace(void);
Point3D RT_tace(Ray ray, double ks, int depth);	//return r, g, b color
void Intersection(bool &h, Ray &r, Point3D &h_color, Point3D &n, Triangle **hit_obj);
void TestShadowEffect(bool &s, Ray &r);