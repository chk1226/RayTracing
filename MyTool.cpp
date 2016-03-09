#include "MyTool.h"


/*******************************************************/
//Triangle

void Triangle::CalculateNormal(void)
{
	Point3D s0, s1;
	s0 = p0 - p1;
	s1 = p2 - p1;

	CrossProduct(n, s0, s1);

}

void Triangle::SetPoints(Point3D p0, Point3D p1, Point3D p2)
{

	this->p0 = p0, this->p1 = p1, this->p2 = p2;
}

void Triangle::Calculate_d(void)
{
	Point3D r = n * p0;
	d = -(r.x+r.y+r.z); 
}

double Triangle::f_of_Q(Point3D q)
{
	return q.x*n.x + q.y*n.y + q.z*n.z + d;
}

bool Triangle::ContainmentTest(Point3D p)
{
	double m0 = CalcualtionLength(p0, p1), 
			m1 = CalcualtionLength(p1, p2), 
			m2 = CalcualtionLength(p0, p2);

	double s = (m0+m1+m2)/2;
	double t_a = sqrt(s*(s-m0)*(s-m1)*(s-m2));	//海龍公式

	double a = CalcualtionLength(p0, p), b = CalcualtionLength(p2, p), c = CalcualtionLength(p1, p);
	s = (a+b+m2)/2;
	double a0 = sqrt(s*(s-a)*(s-b)*(s-m2));
	s = (b+c+m1)/2;
	double a1 = sqrt(s*(s-b)*(s-c)*(s-m1));
	s = (a+c+m0)/2;
	double a2 = sqrt(s*(s-a)*(s-c)*(s-m0));

	if(( (a0+a1+a2)-t_a) > pow(0.1, 6)) return 0;	//大於誤差值 0.1^6
	return 1;

}

void Triangle::SetColor(double r, double g, double b)
{
	color.x = r;
	color.y = g;
	color.z = b;

}

/*******************************************************/
	
/*******************************************************/
//Triangle

void Polygon::SetPoints(const Triangle &t1, const Triangle &t2)
{
	t[0] = t1;
	t[1] = t2;

	t[0].CalculateNormal();
	t[1].CalculateNormal();
	t[0].Calculate_d();
	t[1].Calculate_d();
}



/*******************************************************/




/*******************************************************/
//Other Function

Ray PrimeRay(const Point3D e, const Point3D p)
{
	Ray r;
	r.p = p;
	r.v = p - e;
	r.v = -1.0 * r.v;

	return r;
}

Ray ShadowRay(const Point3D p, const Point3D l)
{
	Ray r;
	r.p = p;
	r.v = l;

	return r;
}

Ray ReflectionRay(const Ray ray, const Point3D n, const Point3D p)
{
	Ray r;
	r.p = p;
	r.v = ray.v - 2*InnerProduct(n,ray.v)*n;
	
	return r;
}

Ray TransmitRay(const Ray ray, const Point3D n, const Point3D p, double kt)
{
	Ray r;
	Point3D L = -1 * ray.v;
	double nl = InnerProduct(n, L);

	r.p = p;
	r.v = ( (1/kt)*nl- sqrt( 1-( 1-pow(nl, 2)/pow(kt, 2) ) ) )*n - (1/kt) * L;

	return r;
}



/*******************************************************/