#include "Render.h"

extern GLbyte *buffer;
extern vector<Obj *> scene_obj;

void DisplayFunc(void)
{
/*	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	
	glBegin(GL_TRIANGLES);
	glColor3f(1,0,0);
	glVertex3f(0.0, 0.0, 40.0);
	glVertex3f(0.0, 40.0, 0.0);	
	glVertex3f(40, 0.0, 0.0);

	glColor3f(0,1,0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 40.0, 0.0);	
	glVertex3f(40.0, 0.0, 0.0);

	glColor3f(0,0,1);
	glVertex3f(0.0, 0.0, 40.0);
	glVertex3f(0.0, 0.0, 0.0);	
	glVertex3f(40.0, 0.0, 0.0);
		
	glColor3f(0,0,0);
	glVertex3f(0.0, 0.0, 40.0);
	glVertex3f(0.0, 40.0, 0.0);	
	glVertex3f(0.0, 0.0, 0.0);
	glEnd();

	glBegin(GL_QUADS);
		glColor3f(0,1,1);

		glVertex3f(40.0, -25.0, 40.0);
		glVertex3f(40.0, -25.0, -40.0);	
		glVertex3f(-40.0, -25.0, -40.0);
		glVertex3f(-40.0, -25.0, 40.0);

	glEnd();
	glBegin(GL_QUADS);
		glColor3f(0,0,1);

		glVertex3f(40.0, -25.0, -40.0);
		glVertex3f(-40.0, -25.0, -40.0);	
		glVertex3f(-40.0, 50.0, -40.0);
		glVertex3f(40.0, 50.0, -40.0);

	glEnd();
	
	glBegin(GL_QUADS);
		glColor3f(0,0,1);

		glVertex3f(20.0, 25.0, 50.0);
		glVertex3f(-20.0, 25.0, 50.0);	
		glVertex3f(-20.0, 50.0, 30.0);
		glVertex3f(20.0, 50.0, 30.0);

	glEnd();

	
	glBegin(GL_LINES);
	glColor3f(1,0,0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(50.0, 0.0, 0.0);
	glColor3f(0,1,0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 50.0, 0.0);
	glColor3f(0,0,1);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 50.0);
	glEnd();
	
*/
	RayTrace();
	cout<<"Render end!!\n";

	glutSwapBuffers();
}

extern const int WIDTH, HEIGHT;
extern const double f_n, f_f, f_l, f_r, f_b, f_t;
extern Point3D e;	//eye postion
extern Point3D u, v, w;

extern const double K_s;

void RayTrace(void)
{
	Point3D c = e - f_n*w;
	Point3D p0 = c + f_l*u + f_b*v,
			p1 = c + f_r*u + f_b*v,
			p2 = c + f_r*u + f_t*v,
			p3 = c + f_l*u + f_t*v;

	Point3D dt_y = (1.0/(HEIGHT-1)) * (p3-p0);
	Point3D dt_x = (1.0/(WIDTH-1)) * (p1-p0);

	Point3D s_p = p0, color;
	Ray r;
	for(int y=0; y<HEIGHT; y++)
	{
		for(int x=0; x<WIDTH; x++)
		{
			s_p = dt_x + s_p;
			r = PrimeRay(e, s_p);
			color = RT_tace(r, K_s, 1);

			int index = ((HEIGHT-1)-y)*(WIDTH*3) + x*3;
			buffer[index] = static_cast<int>(255.0 * color.x);
			buffer[index+1] =  static_cast<int>(255.0 * color.y);
			buffer[index+2] =  static_cast<int>(255.0 * color.z);
		}

		s_p = (static_cast<double>(y) * dt_y) + p0;
	}

	glDrawBuffer(GL_FRONT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 3);			
	glRasterPos2i(0, HEIGHT);
	glDrawPixels(WIDTH, HEIGHT, GL_RGB,GL_UNSIGNED_BYTE, buffer);

}

const int MAX_LEVEL = 7;
const double MIN_KS = 0.2;

extern const Point3D l_pos;
extern const double K_a;
extern const double K_d;

Point3D RT_tace(Ray ray, double ks, int depth)
{
	Point3D z_color = {0.0, 0.0, 0.0},	//zero color
			b_color = {0.0, 0.0, 0.0},	//back ground color
			h_color,					//hit color
			m_color,					//material color
			s_color = {1.0, 1.0, 1.0};	//白色的反光
	Point3D n;
	Triangle *hit_obj = 0;
	bool hit = 0;	//check hit obj

	//check depth
	if(depth>MAX_LEVEL || ks<MIN_KS) return z_color;

	//computer the intersection
	Intersection(hit, ray, m_color, n, &hit_obj);

	if(hit == 0) return b_color;
	h_color = K_a * m_color;

	//check if p is directly lighted
	Point3D p = ray.t*ray.v + ray.p;
	Point3D l = l_pos - p;
	NormalizeV(l);

	if(InnerProduct(n, l) > 0.0)
	{
		Ray s_ray = ShadowRay(p, l);
		bool s = 0;	//是否shadow_ray跟場景發生碰撞, 有回傳1, 否則回傳0
		TestShadowEffect(s, s_ray);
		if(s == 0)
		{
			Point3D v = ::e - p;

			NormalizeV(v);

			h_color = h_color + DiffuseColor(m_color, l, n) + SpecularColor(s_color, v, l, n);

			if(h_color.x>1.0) h_color.x = 1.0;
			if(h_color.y>1.0) h_color.y = 1.0;
			if(h_color.z>1.0) h_color.z = 1.0;
			
		}

	}

	//compute reflection effect
	if(hit_obj->K_s > 0)
	{
		Ray r_ray = ReflectionRay(ray, n, p);
		Point3D r_color = RT_tace(r_ray, K_s*hit_obj->K_s, depth+1);
		h_color = h_color + r_color;

		if(h_color.x>1.0) h_color.x = 1.0;
		if(h_color.y>1.0) h_color.y = 1.0;
		if(h_color.z>1.0) h_color.z = 1.0;
	}

	//transmission effect
	if(hit_obj->K_t>0)
	{
		Ray t_ray = TransmitRay(ray, n, p, hit_obj->K_t);
		Point3D t_color = RT_tace(t_ray, K_s*hit_obj->K_t, depth+1);
		h_color = h_color + t_color;

		if(h_color.x>1.0) h_color.x = 1.0;
		if(h_color.y>1.0) h_color.y = 1.0;
		if(h_color.z>1.0) h_color.z = 1.0;
	}

	return h_color;
}

void Intersection(bool &h, Ray &r, Point3D &h_color, Point3D &n, Triangle **hit_obj)
{
	Triangle *t = 0;
	int index = scene_obj.size();
	double T, reg_T;
	bool f = 1;


	for(int i=0; i<index; i++)
	{
		t = dynamic_cast<Triangle *>(scene_obj[i]);
		if(InnerProduct(t->n, r.v) == 0.0)
		{
			continue;
		}

		T = -(t->f_of_Q(r.p) / InnerProduct(t->n, r.v));


		if(T>pow(0.1,6))	//誤差
		{
			Point3D p = T*r.v + r.p;
						
			if(t->ContainmentTest(p))
			{
				if(f) f = 0, reg_T = T;

				if(T<=reg_T)
				{
					h = 1;
					reg_T = r.t = T;
					n = t->n;
					h_color = t->color;
					*(hit_obj) = t;
				}
			}

		}

	}

}

void TestShadowEffect(bool &s, Ray &r)
{
	Triangle *t = 0;
	int index = scene_obj.size();
	double T;


	for(int i=0; i<index; i++)
	{
		t = dynamic_cast<Triangle *>(scene_obj[i]);
		if(InnerProduct(t->n, r.v) == 0.0)
		{
			continue;
		}

		T = -(t->f_of_Q(r.p) / InnerProduct(t->n, r.v));


		if(T>pow(0.1,6))	//誤差
		{
			Point3D p = T*r.v + r.p;
						
			if(t->ContainmentTest(p))
			{
				s = 1;
				return;
			}

		}

	}

}
