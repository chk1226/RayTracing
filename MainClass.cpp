#include "Render.h"

/*******************************************************/
//Constant

extern const int WIDTH = 800;
extern const int HEIGHT = 800;
extern const double f_n = 1, f_f = 800, f_l = -tan(30.0*ANGLE), f_r = tan(30.0*ANGLE);
extern const double f_b = f_l, f_t = f_r;
Point3D e;	//eye postion
Point3D u, v, w;

GLbyte *buffer = 0;
vector<Obj *> scene_obj;
/*******************************************************/

/*******************************************************/
//Light Parameter
extern const Point3D l_pos = {0.0, 100.0, 0.0};
extern const double K_a = 0.3;
extern const double K_d = 0.6;
extern const double K_s = 0.7;


/*******************************************************/



int offset_x = 120;
int offset_y = 30;

void InitSetting(void)
{

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glOrtho(0.0, WIDTH, HEIGHT, 0.0, 0.0,1.2);
	
	//glEnable(GL_DEPTH_TEST);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();	
	//glFrustum(f_l, f_r, f_b, f_t, f_n, f_f);

	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();

	double r = 200;
	e.x = r*cos(offset_x*ANGLE)*sin(offset_y*ANGLE);
	e.y = r*cos(offset_y*ANGLE);
	e.z = r*sin(offset_y*ANGLE)*sin(offset_x*ANGLE);

	//glViewport(0, 0, WIDTH, HEIGHT);
	//gluLookAt(  e.x, e.y, e.z,
	//			0.0, 0.0, 0.0,
	//			0.0, 1.0, 0.0);
	//
	//glClearColor(1.0, 1.0, 1.0, 1.0);

	w = 0.0 - e; NormalizeV(w); 
	v.x = v.z = 0.0, v.y = 1.0;
	CrossProduct(u, v, w);
	CrossProduct(v, w, u);

}

void DeleteAllObject(void);
void Keyboard(unsigned char c, int x, int y)
{
	switch(c)
	{
		case 0x1b: 
			DeleteAllObject();
			exit(0);
		case 'a':
		case 'A': 
			if(offset_x<350) offset_x += 10;
			break;
		case 'd':
		case 'D': 
			if(offset_x>10) offset_x -= 10;
			break;
		case 's':
		case 'S': 
			if(offset_y<170) offset_y += 10;
			break;
		case 'w':
		case 'W': 
			if(offset_y>10) offset_y -= 10;
			break;

	}
	InitSetting();
	glutPostRedisplay();
}

void InitObject(void)
{
	buffer = new GLbyte[WIDTH*HEIGHT*3];

	//三角錐
	Point3D p0 = {0.0, 0.0, 40.0},
			p1 = {0.0, 40.0, 0.0},
			p2 = {40.0, 0.0, 0.0};
	Triangle *ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(1.0, 0, 0);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.0;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);

	p0.z = 0.0;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0, 1.0, 0);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.0;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);

	p0.z = 40.0;
	p1.y = 0.0;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0, 0, 1.0);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.0;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);

	p1.y = 40.0;
	p2.x = 0.0;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(1.0, 1.0, 0);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.0;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);

	//地板
	p0.x = -40, p0.y = -25, p0.z = -40;
	p1.x = 40, p1.y = -25, p1.z = -40;
	p2.x = 40, p2.y = -25, p2.z = 40;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0.2, 0.2, 0.2);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.5;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);
	p0.x = 40, p0.y = -25, p0.z = 40;
	p1.x = -40, p1.y = -25, p1.z = 40;
	p2.x = -40, p2.y = -25, p2.z = -40;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0.2, 0.2, 0.2);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.5;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);

	p0.x = 40, p0.y = -25, p0.z = -41;
	p1.x = -40, p1.y = -25, p1.z = -41;
	p2.x = -40, p2.y = 50, p2.z = -41;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0.2, 0.2, 0.2);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.5;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);
	p0.x = -40, p0.y = 50, p0.z = -41;
	p1.x = 40, p1.y = 50, p1.z = -41;
	p2.x = 40, p2.y = -25, p2.z = -41;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0.2, 0.2, 0.2);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.5;
	ptr->K_t = 0.0;
	scene_obj.push_back(ptr);

	//透明遮蔽物
	p0.x = 20, p0.y = 25, p0.z = 50;
	p1.x = -20, p1.y = 25, p1.z = 50;
	p2.x = -20, p2.y = 50, p2.z = 30;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0.4, 0.4, 0.4);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.5;
	ptr->K_t = 0.5;
	scene_obj.push_back(ptr);
	p0.x = -20, p0.y = 50, p0.z = 30;
	p1.x = 20, p1.y = 50, p1.z = 30;
	p2.x = 20, p2.y = 25, p2.z = 50;
	ptr = new Triangle;
	ptr->SetPoints(p0, p1, p2);
	ptr->SetColor(0.4, 0.4, 0.4);
	ptr->CalculateNormal();
	ptr->Calculate_d();
	ptr->K_s = 0.5;
	ptr->K_t = 0.5;
	scene_obj.push_back(ptr);

}

void DeleteAllObject(void)
{
	delete [] buffer;

	int index = scene_obj.size();
	for(int i=0; i<index; i++)
	{
		delete scene_obj[i];
	}
	
	scene_obj.clear();
}

int main(int argc, char *argv[])
{

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutCreateWindow("Ray-Tracing");
	InitObject();
	InitSetting();
	glutDisplayFunc(DisplayFunc);
	glutKeyboardFunc(Keyboard);

	glutMainLoop();
	return 0;
}
