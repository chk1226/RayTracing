#include "MyMath.h"

Point3D operator-(const Point3D &a, const Point3D &b)
{
	Point3D reg = {a.x-b.x, a.y-b.y, a.z-b.z};

	return reg;
}

Point3D operator-(const double v, const Point3D &b)
{
	Point3D reg = {v-b.x, v-b.y, v-b.z};

	return reg;
}

Point3D operator*(const double a, const Point3D &b)
{
	Point3D reg = {a*b.x, a*b.y, a*b.z};
	
	return reg;
}

Point3D operator*(const Point3D &a, const Point3D &b)
{
	Point3D reg = {a.x*b.x, a.y*b.y, a.z*b.z};

	return reg;
}

Point3D operator+(const Point3D &a, const Point3D &b)
{
	Point3D reg = {a.x+b.x, a.y+b.y, a.z+b.z};

	return reg;
}

void NormalizeV(Point3D &t)
{
	double length;
	
	length = sqrt(pow(t.x,2) + pow(t.y,2) + pow(t.z, 2));

	if(length == 0.0f)
		length = 1.0f;

	t.x /= length, t.y /= length ,t.z /= length;

}

double InnerProduct(const Point3D &v1, const Point3D &v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

void CrossProduct(Point3D &u, const Point3D &w, const Point3D &view_up)
{
	u.x = w.y*view_up.z - w.z*view_up.y;
	u.y = w.z*view_up.x - w.x*view_up.z;
	u.z = w.x*view_up.y - w.y*view_up.x;

	NormalizeV(u);
}

double CalcualtionLength(const Point3D &p0, const Point3D &p1)
{
	return sqrt( pow(p0.x-p1.x,2)+pow(p0.y-p1.y,2)+pow(p0.z-p1.z,2) ) ;
}

extern const double K_d;
extern const double K_s;
Point3D DiffuseColor(Point3D color, const Point3D l, const Point3D n)
{
	double d = 0.0;
	if( (d = InnerProduct(n, l)) < 0) d = 0.0;

	color = K_d*d*color;

	return color;
}

Point3D SpecularColor(Point3D color, const Point3D v, const Point3D l, const Point3D n)
{

	Point3D r;
	double d = 2*InnerProduct(n, l);

	r = d * n - l;
	NormalizeV(r);

	d = InnerProduct(r, v); 
	if(d<0) d = 0.0;

	d = pow(d, 32);

	color = K_s*d*color;

	return color;
}



/*
void Gauss(double *value, double *cons, int dimension)
{
	double c = 0.0;
	for(int j=0; j<dimension; j++)
	{
		for(int i=0;i<dimension-1;i++)
			for(int k=j;k<dimension-1;k++)
			{
				c = value[k*dimension+j];	
				if(c == 0)
				{
					for(int l=0; l<dimension; l++)
						Swap(value[k*dimension + l], value[(k+1)*dimension + l]);
					Swap(cons[k], cons[1+k]);
				}
			}

		c = value[j*dimension+j];	

		for(int i=0; i<dimension; i++)
		{
			value[j*dimension+i]/=c;
		}
		cons[j] /= c;
		
		for(int i=0; i<dimension; i++)
		{
			if(i!=j)
			{
				c = value[i*dimension+j];
				for(int k=0; k<dimension; k++)
				{
					value[i*dimension+k]-= c*value[j*dimension+k];
				}
				cons[i]-= c*cons[j];
			}
		}
	}

}
*/