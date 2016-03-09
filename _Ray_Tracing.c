#include<stdio.h>
#include<math.h>
#include<GL/glut.h>
#include<stdlib.h>
#include<math.h>


#define WSIZE 800
#define MAX_LEVEL 5
#define MIN_KS 0.2
#define TRUE 0
#define FALSE 1
#define POLY 1
#define SPHERE 0
#define PI 3.141592653

unsigned char color_buffer[WSIZE][WSIZE][3],Z_buffer[WSIZE][WSIZE][3];
double normal_buffer[WSIZE][WSIZE][3];


struct GeoMat{
	struct GeoMat *last;
	double Matrix[4][4];
}GeoStack,NowMat; 
double ViewMat[4][4],ProMat[4][4],ScrMat[2][2];

struct Poly{
	double vertex[3][3];
	unsigned char color[3];
	double normal[3];
	double unfaithful[3];//み
	double radius;
	double Ks,Kt,reIndex;
}**PolyNull;

struct Sphere{
	double center[3];
	unsigned char color[3];
	double radius;
	double Ks,Kt,reIndex;
}*SphereNull;

struct Model{
	int Type;
	int count;
	struct Poly **subset;
	struct Sphere *S;
}*ModelNull;

struct Object{
	int index;
	struct Object *next;
	struct Model *model;
}*scene,*ObjectNull;

struct Light{
	double lightColor[3];
	double lightSet[3];
	struct Light *next;
}*lightSourceHead,*lightNull;


double  points[][3] = {{-1, -1, -1}, {1, -1, -1}, {1,-1,1}, {-1,-1, 1}, 
						{-1,1, -1}, {1,1, -1},{1,1,1}, {-1,1,1}};
int  face[][4] = {{0, 1, 2, 3}, {7,6,5,4}, {0,4,5,1}, 
                    {1, 5, 6, 2}, {2, 6, 7 , 3}, {3, 7, 4, 0}};
/*view value*/
double Center[3] ={0,0,0};
double eye[3] = {4,7,6},EyeRota = 0;
double l,r,t,b,n,f;
/*---light value---*/
double Ka=0.3,Kd=0.9,Ks=0.7;
int shininess = 32;
double lightSourceCount=0,LS[3]={1,8,1};
unsigned char ambient_color[3] = {255,255,255},diffuse_color[3] = {255,255,255},speculor_color[3] ={255,255,255};



void PushM()
{
	struct GeoMat *temp;
	temp = (struct GeoMat *)malloc(sizeof(struct GeoMat));
	*temp = NowMat;
	NowMat.last = temp;

}
void PopM()
{
	Center[0] = 0;
	Center[1] = 0;
	Center[2] = 0;
	NowMat = *NowMat.last;
}
void printfM(double M[4][4])
{
	int i,j;
	for(i=0;i<4;i++){
		printf("|");
		for(j=0;j<4;j++){
			printf(" %f",M[i][j]);
		}
		printf(" |\n");
	}
	printf("\n");
	
}
void productM(double M1[4][4],double M2[4][4],double result[4][4])
{
	int i,j;
	double temp[4][4];
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			temp[i][j] = M1[i][0]*M2[0][j]+M1[i][1]*M2[1][j]+M1[i][2]*M2[2][j]+M1[i][3]*M2[3][j];
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			result[i][j] = temp[i][j];
}
double* productM_P(double x,double y,double z,int v)
{
	double *result;
	int i;
	result = (double *)malloc(sizeof(double)*3);
	for(i=0;i<4;i++)
		result[i] = x*NowMat.Matrix[i][0] + y*NowMat.Matrix[i][1] +
		            z*NowMat.Matrix[i][2] + v*NowMat.Matrix[i][3];
	return result;
}
void set_light_source(double x,double y,double z,double r,double g,double b){
	struct Light *light,*now;
	now = lightSourceHead;
	while(now->next != lightNull )
		now = now->next;
	light = (struct Light *)malloc(sizeof(struct Light));
	light->lightColor[0] = r;
	light->lightColor[1] = g;
	light->lightColor[2] = b;
	light->lightSet[0] = x;
	light->lightSet[1] = y;
	light->lightSet[2] = z;
	light->next = lightNull;

	now->next =light;
	lightSourceCount ++;
}

void init_light(){

	double *lightSet;
	lightSourceHead = (struct Light *)malloc(sizeof(struct Light));
	lightNull = (struct Light *)malloc(sizeof(struct Light));
	lightSourceHead->next = lightNull;
	lightSourceCount = 0;

//	Rotate(LightAngle,0,1,0);
	
	lightSet = productM_P(LS[0],LS[1],LS[2],1);
	set_light_source(lightSet[0],lightSet[1],lightSet[2],1,1,1);

}



void init_object_struct()
{
	PolyNull = (struct Poly **)malloc(sizeof(struct Poly));

	ModelNull = (struct Model *)malloc(sizeof(struct Model));
	ModelNull->count = 0;
	ModelNull->subset = PolyNull;
	
	ObjectNull = (struct Object *)malloc(sizeof(struct Object));
	ObjectNull->index = -1;
	ObjectNull->model = ModelNull;
	ObjectNull->next = ObjectNull;

	scene = (struct Object *)malloc(sizeof(struct Object));
	scene->index = 0;
	scene->model = ModelNull;
	scene->next = ObjectNull;
}
void init_func()
{
	int i,j,c;
	glClearColor(0.0, 0.0, 0.0, 1.0);      /*set the background color BLACK */
                   
	glClear(GL_COLOR_BUFFER_BIT);		/*Clear the Depth & Color Buffers */
	for(i=0;i<WSIZE;i++)
		for(j=0;j<WSIZE;j++)
			for(c=0;c<3;c++){
				color_buffer[i][j][c] = 0;	
			}
	for(i=0;i<WSIZE;i++)
		for(j=0;j<WSIZE;j++)
			for(c=0;c<3;c++){
				normal_buffer[i][j][c] = 0;	
			}
	/*init Geometrical Matrix*/
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			if(i==j)
				GeoStack.Matrix[i][j] = 1;
			else
				GeoStack.Matrix[i][j] = 0;
	GeoStack.last = NULL;
	NowMat = GeoStack;
	init_object_struct();
	
	/*init light source head*/
	PushM();
	init_light();
	PopM();

	glFlush();
}

void normalize(double x[3],double result[3]){	//vector normalize
	double temp;
	temp = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	if(temp!=0){
		result[0] = x[0]/temp;
		result[1] = x[1]/temp;
		result[2] = x[2]/temp;
	}else result[0] = result[1] = result[2] = 0;
}

void cross(double A[3],double Bx,double By,double Bz,double C[3])	//vector A Cross vector B = C vector
{
	C[0] = A[1]*Bz - A[2]*By;
	C[1] = A[2]*Bx - A[0]*Bz;
	C[2] = A[0]*By - A[1]*Bx;
}


void set_poly_normal_unfaithful(struct Poly *P)
{
	int i;
	double ABvector[3],ACvector[3],normal[3];
	double funcA[3],funcB[3],funcC[3],funcD[3];//0:ABキだよ祘Α把计 1:ACキだよ祘Α把计 2:ABCキよ祘Α把计
	double ABmedian[3],ACmedian[3];
	double det,Upoint[3];
	for(i=0;i<3;i++){
		ABvector[i] = P->vertex[1][i] - P->vertex[0][i];
		ACvector[i] = P->vertex[2][i] - P->vertex[0][i];
	}
	normalize(ABvector,ABvector);
	normalize(ACvector,ACvector);
	cross(ABvector,ACvector[0],ACvector[1],ACvector[2],normal);
	normalize(normal,normal);
	/*---―ABい翴 & ACい翴*/
	for(i=0;i<3;i++){
		ABmedian[i] = (P->vertex[0][i] + P->vertex[1][i])/2;
		ACmedian[i] = (P->vertex[0][i] + P->vertex[2][i])/2;
	}
	/*---―ABキだよ祘Α---*/	
	funcA[0] = ABvector[0];
	funcB[0] = ABvector[1];
	funcC[0] = ABvector[2];
	funcD[0] = funcA[0]*ABmedian[0] + funcB[0]*ABmedian[1] + funcC[0]*ABmedian[2];
	/*---―ACキだよ祘Α---*/
	funcA[1] = ACvector[0];
	funcB[1] = ACvector[1];
	funcC[1] = ACvector[2];
	funcD[1] = funcA[1]*ACmedian[0] + funcB[1]*ACmedian[1] + funcC[1]*ACmedian[2];
	/*---―ABCキよ祘Α---*/
	funcA[2] = normal[0];
	funcB[2] = normal[1];
	funcC[2] = normal[2];
	funcD[2] = funcA[2]*P->vertex[0][0] + funcB[2]*P->vertex[0][1] + funcC[2]*P->vertex[0][2]; 
	/*---―み---*/	
	//[funcA[0] funcB[0] funcC[0]]   [x]   [funcD[0]]
	//[funcA[1] funcB[1] funcC[1]] * [y] = [funcD[1]]
	//[funcA[2] funcB[2] funcC[2]]   [z]   [funcD[2]]
	det = funcA[0]*funcB[1]*funcC[2] + funcB[0]*funcC[1]*funcA[2] + funcC[0]*funcA[1]*funcB[2] 
		- funcC[0]*funcB[1]*funcA[2] - funcB[0]*funcA[1]*funcC[2] - funcA[0]*funcC[1]*funcB[2];

	Upoint[0] = (  (funcB[1]*funcC[2] - funcB[2]*funcC[1])*funcD[0] 
		         - (funcB[0]*funcC[2] - funcB[2]*funcC[0])*funcD[1]
				 + (funcB[0]*funcC[1] - funcB[1]*funcC[0])*funcD[2])/det;
	Upoint[1] = (- (funcA[1]*funcC[2] - funcA[1]*funcC[1])*funcD[0] 
		         + (funcA[0]*funcC[2] - funcA[2]*funcC[0])*funcD[1]
				 - (funcA[0]*funcC[1] - funcA[1]*funcC[0])*funcD[2])/det;
	Upoint[2] = (  (funcA[1]*funcB[2] - funcA[2]*funcB[1])*funcD[0]
				 - (funcA[0]*funcB[2] - funcA[2]*funcB[0])*funcD[1]
				 + (funcA[0]*funcB[1] - funcA[1]*funcB[0])*funcD[2])/det;
	/*---―畖---*/
	P->radius =sqrt( pow(Upoint[0] - P->vertex[0][0],2) + pow(Upoint[1] - P->vertex[0][1],2) + pow(Upoint[2] - P->vertex[0][2],2));
	for(i=0;i<3;i++){
		P->normal[i] = normal[i];
		P->unfaithful[i] = Upoint[i];
	}
	
}

struct Poly *set_poly(double vertex0[3],double vertex1[3],double vertex2[3],unsigned char color[3])
{	// set polygon vector and color and normal
	int xyz,i;
	struct Poly *P;
	
	P = (struct Poly *)malloc(sizeof(struct Poly));
	for(xyz=0;xyz<3;xyz++){
		P->vertex[0][xyz] = vertex0[xyz];
		P->vertex[1][xyz] = vertex1[xyz];
		P->vertex[2][xyz] = vertex2[xyz];
	}
	set_poly_normal_unfaithful(P);

	for(i=0;i<3;i++)
		P->color[i] = color[i];
	
	return P;
}
struct Sphere *set_sphere(double center[3],double radius,unsigned char color[3])
{
	int i;
	struct Sphere *S;
	S = (struct Sphere*)malloc(sizeof(struct Sphere));
	for(i=0;i<3;i++){
		S->center[i] = center[i];
		S->color[i] = color[i];
	}
	S->radius = radius;
	return S;
}
struct Model *set_poly_model(struct Poly **poly,int PolyCount)
{
	struct Model *M;
	
	M = (struct Model *)malloc(sizeof(struct Model));	
	M->Type = POLY;
	M->count = PolyCount;		
	M->subset = poly;
	M->S = SphereNull;
	return M;
}
struct Model *set_sphere_model(struct Sphere *S)
{
	struct Model *M;
	M = (struct Model *)malloc(sizeof(struct Model));	
	M->Type = SPHERE;
	M->count = 1;
	M->subset = PolyNull;
	M->S = S;
	return M;
}
void set_object(struct Model *floorM)
{
	struct Object *O,*Temp;
	O = (struct Object *)malloc(sizeof(struct Object));
	Temp = scene;
	while(Temp->next!=ObjectNull){
		Temp = Temp->next;
	}
	O->model = floorM;
	O->next = ObjectNull;
	O->index = Temp->index + 1;
	Temp->next = O;
}

void draw_wall(double Ks,double Kt,double reIndex,int type)
{
	struct Poly **wallP;
	struct Model *wallM;
	double wallV[4][3] = {{-5,0,0},{5,0,0},{5,5,0},{-5,5,0}};
	double *w[4];
	unsigned char Vcolor[4][2][3] = {{{255,255,255},{255,255,255}},{{0,0,0},{0,0,0}}
									,{{255,0,0},{255,0,0}},{{0,125,0},{0,125,0}}};
	int PolyCount = 2;
	wallP = (struct Poly **)malloc(sizeof(struct Poly)*PolyCount);
	w[0] = productM_P(wallV[0][0],wallV[0][1],wallV[0][2],1);	//four node in word coordinate
	w[1] = productM_P(wallV[1][0],wallV[1][1],wallV[1][2],1);
	w[2] = productM_P(wallV[2][0],wallV[2][1],wallV[2][2],1);
	w[3] = productM_P(wallV[3][0],wallV[3][1],wallV[3][2],1);
	wallP[0] = set_poly(w[0],w[1],w[2],Vcolor[type][0]);
	wallP[1] = set_poly(w[2],w[3],w[0],Vcolor[type][1]);
	wallP[0]->Ks=Ks;
	wallP[1]->Ks=Ks;
	wallP[0]->Kt=Kt;
	wallP[1]->Kt=Kt;
	wallP[0]->reIndex = reIndex;
	wallP[1]->reIndex = reIndex;
	wallM = set_poly_model(wallP,PolyCount);
	set_object(wallM);
}

void draw_floor(double Ks,double Kt,double reIndex)
{
	struct Poly **floorP;
	struct Model *floorM;
	double floorV[4][3] = {{-5,0,5},{5,0,5},{5,0,-5},{-5,0,-5}};
	double *w[4];
	unsigned char Vcolor[2][3] = {{255,255,255},{255,255,255}};
	int PolyCount = 2;
	floorP = (struct Poly **)malloc(sizeof(struct Poly)*PolyCount);
	w[0] = productM_P(floorV[0][0],floorV[0][1],floorV[0][2],1);	//four node in word coordinate
	w[1] = productM_P(floorV[1][0],floorV[1][1],floorV[1][2],1);
	w[2] = productM_P(floorV[2][0],floorV[2][1],floorV[2][2],1);
	w[3] = productM_P(floorV[3][0],floorV[3][1],floorV[3][2],1);

	floorP[0] = set_poly(w[0],w[1],w[2],Vcolor[0]);
	floorP[1] = set_poly(w[2],w[3],w[0],Vcolor[1]);
	floorP[0]->Ks=Ks;
	floorP[1]->Ks=Ks;
	floorP[0]->Kt=Kt;
	floorP[1]->Kt=Kt;
	floorP[0]->reIndex = reIndex;
	floorP[1]->reIndex = reIndex;
	floorM = set_poly_model(floorP,PolyCount);
	set_object(floorM);

}
void draw_cube(double Ks,double Kt,double reIndex)
{
	struct Poly **cubeP;
	struct Model *cubeM;
	double *w[4];
	unsigned char Vcolor[6][3] = {{255,0,0},{0,255,0},{0,0,255},{255,255,0},{255,0,255},{0,255,255}};
	int PolyCount = 12,PolyNum;

	cubeP = (struct Poly **)malloc(sizeof(struct Poly)*PolyCount);
	for(PolyNum = 0 ; PolyNum < PolyCount ; PolyNum += 2){
		w[0] = productM_P(points[face[PolyNum/2][0]][0],points[face[PolyNum/2][0]][1],points[face[PolyNum/2][0]][2],1);	//four node in word coordinate
		w[1] = productM_P(points[face[PolyNum/2][1]][0],points[face[PolyNum/2][1]][1],points[face[PolyNum/2][1]][2],1);
		w[2] = productM_P(points[face[PolyNum/2][2]][0],points[face[PolyNum/2][2]][1],points[face[PolyNum/2][2]][2],1);
		w[3] = productM_P(points[face[PolyNum/2][3]][0],points[face[PolyNum/2][3]][1],points[face[PolyNum/2][3]][2],1);	
		cubeP[PolyNum] = set_poly(w[0],w[1],w[2],Vcolor[PolyNum/2]);
		cubeP[PolyNum+1] = set_poly(w[2],w[3],w[0],Vcolor[PolyNum/2]);
		cubeP[PolyNum]->Ks = Ks;
		cubeP[PolyNum+1]->Ks = Ks;
		cubeP[PolyNum]->Kt = Kt;
		cubeP[PolyNum+1]->Kt = Kt;
		cubeP[PolyNum]->reIndex = reIndex; 
		cubeP[PolyNum+1]->reIndex = reIndex;
	}
	

	cubeM = set_poly_model(cubeP,PolyCount);
	set_object(cubeM);
}
void draw_ball(double Ks,double Kt,double reIndex)
{
	struct Sphere *ballS;
	struct Model *ballM;
	double *w;
	double center[3] = {0,0,0};
	unsigned char Vcolor[3] = {255,255,255};
	
	ballS = (struct Sphere *)malloc(sizeof(struct Sphere));
	w = productM_P(center[0],center[1],center[2],1);	
			
	ballS = set_sphere(w,1,Vcolor);
	ballS->Ks = Ks;
	ballS->Kt = Kt;
	ballS->reIndex = reIndex;
	ballM = set_sphere_model(ballS);
	set_object(ballM);
}
void Scale(double sx,double sy,double sz)
{
	double SM[4][4];
	int i,j;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
				SM[i][j] = 0;
	SM[0][0] = sx;
	SM[1][1] = sy;
	SM[2][2] = sz;
	SM[3][3] = 1;
	
	productM(NowMat.Matrix,SM,NowMat.Matrix);
	
}

void Translate(double dx,double dy,double dz)
{
	double TM[4][4];
	int i,j;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			if(i==j)
				TM[i][j] = 1;
			else
				TM[i][j] = 0;
	TM[0][3] = dx;
	TM[1][3] = dy;
	TM[2][3] = dz;
	Center[0] += dx;
	Center[1] += dy;
	Center[2] += dz;
	productM(TM,NowMat.Matrix,NowMat.Matrix);
}
void T(double dx,double dy,double dz,double TM[4][4])
{
	int i,j;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			if(i==j)
				TM[i][j] = 1;
			else
				TM[i][j] = 0;
	TM[0][3] = dx;
	TM[1][3] = dy;
	TM[2][3] = dz;
}

void Rx(int N,double rx[4][4],double dx,double dy,double dz)
{
	
	int i,j;
	double d;
	
	d = sqrt(dz*dz+dy*dy);
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
				rx[i][j] = 0;
	rx[0][0] = 1;
	if(d!=0){
		rx[1][1] = dz/d;
		rx[1][2] = -dy/d;
		rx[2][1] = dy/d;
		rx[2][2] = dz/d;
	}else{
		rx[1][1] = 1;
		rx[2][2] = 1;
	}
	rx[3][3] = 1;
	if(N == -1){
		rx[1][2] *= -1;
		rx[2][1] *= -1;
	}

} 
void Ry(int N,double ry[4][4],double dx,double dy,double dz)
{
	
	int i,j;
	double d;
	
	d = sqrt(dz*dz+dy*dy);
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
				ry[i][j] = 0;
	ry[0][0] = d;
	ry[0][2] = dx;
	ry[1][1] = 1;
	ry[2][0] = -dx; 
	ry[2][2] = d;
	ry[3][3] = 1;
	if(N == -1){
		ry[0][2] *= -1;
		ry[2][0] *= -1;
	}
} 
void Rz(double drgee,double rz[4][4])
{
	
	int i,j;

	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
				rz[i][j] = 0;
	rz[0][0] = cos(drgee*PI/180.0);
	rz[0][1] = -sin(drgee*PI/180.0);
	rz[1][0] = sin(drgee*PI/180.0); 
	rz[1][1] = cos(drgee*PI/180.0);
	rz[2][2] = 1;
	rz[3][3] = 1;

}

void Rotate(double drgee,double x,double y,double z)
{
	double rx[4][4],ry[4][4],rz[4][4],TM[4][4];

	
	T(-Center[0],-Center[1],-Center[2],TM);//T(-P)
	productM(TM,NowMat.Matrix,NowMat.Matrix);
	
	Rx(1,rx,x,y,z);//Rx(sitaX)
	productM(rx,NowMat.Matrix,NowMat.Matrix);
	Ry(-1,ry,x,y,z);//Ry(-sitaY)
	productM(ry,NowMat.Matrix,NowMat.Matrix);
	Rz(drgee,rz);//Rz(sita)
	productM(rz,NowMat.Matrix,NowMat.Matrix);
	Ry(1,ry,x,y,z);//Ry(sitaY)
	productM(ry,NowMat.Matrix,NowMat.Matrix);
	Rx(-1,rx,x,y,z);//Rx(-sitaX)
	productM(rx,NowMat.Matrix,NowMat.Matrix);


	T(Center[0],Center[1],Center[2],TM);//T(P)	
	productM(TM,NowMat.Matrix,NowMat.Matrix);
	
//	printfM(NowMat.Matrix);
}
void draw()
{

	PushM();
	Translate(0,0,-5);
	draw_wall(0,0,1,0);
	PopM();
	PushM();
	Translate(-5,0,0);
	Rotate(90,0,1,0);
	draw_wall(1,0,1,1);
	PopM();
	PushM();
	Translate(5,0,0);
	Rotate(270,0,1,0);
	draw_wall(0,0,1,2);
	PopM();
	PushM();
	draw_floor(0,0,1);
	PopM();
	PushM();
	Translate(-1,2,-1);
	Scale(2,2,2);
	draw_cube(0.2,0,1);
	PopM();
	PushM();
	Translate(-1,1,2);
	draw_ball(0.5,0,1);
	PopM();
	PushM();
	Translate(0,0,5);
	draw_wall(0,0.9,1,3);
	PopM();
}





void LookAt(float Ex,float Ey,float Ez,float Fx,float Fy,float Fz,float Vx,float Vy,float Vz)
{
	double u[3],v[3],nw[3],M[4][4],t[4][4],fabsU,fabsV,fabsNW;
	int i ,j;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			M[i][j] = 0;

	
	nw[0] = (Fx-Ex);
	nw[1] = (Fy-Ey);
	nw[2] = (Fz-Ez);
	fabsNW = fabs(sqrt(nw[0]*nw[0] + nw[1]*nw[1] + nw[2]*nw[2]));
	if(fabsNW != 0){
		nw[0] = nw[0]/fabsNW;
		nw[1] = nw[1]/fabsNW;
		nw[2] = nw[2]/fabsNW;
	}else{
		nw[0] = nw[1] = nw[2] =0;
	}
	cross(nw,Vx,Vy,Vz,u);
	fabsU = fabs(sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
	if(fabsU!=0){
		u[0] = u[0]/fabsU;
		u[1] = u[1]/fabsU;
		u[2] = u[2]/fabsU;
	}else{
		u[0] = u[1] = u[2] = 0;
	}

	cross(u,nw[0],nw[1],nw[2],v);
	fabsV = fabs(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
	if(fabsV!=0){
	v[0] = v[0]/fabsV;
	v[1] = v[1]/fabsV;
	v[2] = v[2]/fabsV;
	}else{
		v[0] = v[1] = v[2] = 0;
	}
	for(i = 0;i<3;i++){
		M[0][i] = u[i];
		M[1][i] = v[i];
		M[2][i] = -nw[i];
	}
	M[3][3] = 1;
	T(-Ex,-Ey,-Ez,t);
	productM(M,t,ViewMat);
	productM(NowMat.Matrix,ViewMat,ViewMat);
	printf("Viewing Matrix:\n");
	printfM(ViewMat);
	
}
void projection(float xmin,float ymin,float zmin,float xmax,float ymax,float zmax)
{
	int i,j;
//	float l,r,b,t,n,f;
	l = xmin;
	r = xmax;
	b = ymin;
	t = ymax;
	n = zmin;
	f = zmax;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			ProMat[i][j] = 0;
	ProMat[0][0] = 2*n/(r-l);
	ProMat[0][2] = (r+l)/(r-l);
	ProMat[1][1] = 2*n/(t-b);
	ProMat[1][2] = (t+b)/(t-b);
	ProMat[2][2] = -(f+n)/(f-n);
	ProMat[2][3] = (-2)*f*n/(f-n);
	ProMat[3][2] = -1;
	printf("Projection Matrix:\n");
	printfM(ProMat);
}



double triangle_area(double point0[2],double point1[2],double point2[2])
{
	double result;
	result = (point0[0]*point1[1]+point1[0]*point2[1]+point2[0]*point0[1]
		     -point1[0]*point0[1]-point2[0]*point1[1]-point0[0]*point2[1])/2;
	return result;
}

int intersect_sphere(struct Sphere *TestSphere,double Q[3],double d[3],double HitPoint[3])
{
	double A,B,C,D;
	double t[2];
	int n;
	A = 1;
	B = 2*(d[0]*(Q[0] - TestSphere->center[0]) + d[1]*(Q[1] - TestSphere->center[1]) + d[2]*(Q[2] -TestSphere->center[2]));
	C = pow(Q[0] - TestSphere->center[0],2) + pow(Q[1] - TestSphere->center[1],2) + pow(Q[2] - TestSphere->center[2],2) - pow(TestSphere->radius,2);
	D = pow(B,2) - 4*A*C;
	
	if(D <= 0)
		return FALSE;
	else{
		t[0] = (-B+sqrt(D))/2;
		t[1] = (-B-sqrt(D))/2;
		
		if(t[0]<t[1]){
			if(t[0]>0.01)
				n=0;
			else if(t[1]>0.01)
				n=1;
			else
				return FALSE;
		}else if(t[0]>t[1]){
			if(t[1]>0.01)
				n=1;
			else if(t[0]>0.01)
				n=0;
			else
				return FALSE;
		}
		
		HitPoint[0] = Q[0] + d[0]*t[n];
		HitPoint[1] = Q[1] + d[1]*t[n];
		HitPoint[2] = Q[2] + d[2]*t[n];
		return TRUE;
	}
}

int intersect_poly(struct Poly *TestPoly,double Q[3],double d[3],double HitPoint[3])
{
	float FuncA,FuncB,FuncC,FuncD,t,N_d;
//	float A,B,C,D;
	int i,U,V;
	double vertexUV[3][2],Ruv[2];
	double sTra[3],bTra,temp;

	FuncA = TestPoly->normal[0];
	FuncB = TestPoly->normal[1];
	FuncC = TestPoly->normal[2];
	FuncD = -1*(FuncA * TestPoly->vertex[0][0] + FuncB * TestPoly->vertex[0][1] + FuncC * TestPoly->vertex[0][2]);
	N_d = FuncA*d[0] + FuncB*d[1] + FuncC*d[2];//normal dot d vector
	
	if(N_d >= 0)	// no intersection or inisible
		return FALSE;
	t = -1*(FuncA*Q[0] + FuncB*Q[1] + FuncC*Q[2] + FuncD)/(N_d);
	if(t <= 0)
		return FALSE;

	/*---containment test---*/
//	x = Q[0] + d[0]*t;
//	y = Q[1] + d[1]*t;
//	z = Q[2] + d[2]*t;

	
	if(fabs(TestPoly->normal[0]) > fabs(TestPoly->normal[1]))
		if(fabs(TestPoly->normal[0]) > fabs(TestPoly->normal[2])){
			U = 1;
			V = 2;
		}else{
			U = 0;
			V = 1;
		}
	else
		if(fabs(TestPoly->normal[1]) > fabs(TestPoly->normal[2])){
			U = 0;
			V = 2;
		}else{
			U = 0;
			V = 1;
		}
	
	Ruv[0] = Q[U] + d[U]*t;
	Ruv[1] = Q[V] + d[V]*t;
	for(i=0;i<3;i++){
		vertexUV[i][0] = TestPoly->vertex[i][U];
		vertexUV[i][1] = TestPoly->vertex[i][V];
	}
	bTra = fabs(triangle_area(vertexUV[0],vertexUV[1],vertexUV[2]));
	sTra[0] = fabs(triangle_area(Ruv,vertexUV[1],vertexUV[2]));
	sTra[1] = fabs(triangle_area(vertexUV[0],Ruv,vertexUV[2]));
	sTra[2] = fabs(triangle_area(vertexUV[0],vertexUV[1],Ruv));
	temp =fabs( bTra - sTra[0] - sTra[1] - sTra[2])/bTra;
	
	if(temp <= 0.0001 && t > 0.001){
		HitPoint[0] = Q[0] + d[0]*t;
		HitPoint[1] = Q[1] + d[1]*t;
		HitPoint[2] = Q[2] + d[2]*t;
		return TRUE;
	}else
		return FALSE;
		
}

int get_hit(double d[3],double Q[3],struct Object *Obj,struct Poly *HitFace,double HitPoint[3])
{
	int i,result = FALSE;
	struct Object *ObjTemp;
	struct Model *ModTemp;
	struct Poly *PolyTemp;
	struct Sphere *SphereTemp;
	double nearDistance,tempDistance,N[3];
	PolyTemp = (struct Poly *)malloc(sizeof(struct Poly));
	ObjTemp = Obj;
	nearDistance = f;
	while(Obj->next != ObjectNull){
		Obj = Obj->next;
		ModTemp = Obj->model;
		if(ModTemp != ModelNull){
			if(ModTemp->Type == POLY){
				if(ModTemp->subset != PolyNull){
					for(i = 0; i < ModTemp->count; i++){
						PolyTemp = ModTemp->subset[i];
												
						if(intersect_poly(PolyTemp,Q,d,HitPoint) == TRUE){
							tempDistance = sqrt(pow(Q[0] - HitPoint[0],2) + pow(Q[1] - HitPoint[1],2) + pow(Q[2] - HitPoint[2],2));
							if(nearDistance > tempDistance){
								*HitFace = *PolyTemp;
								result = TRUE;
								nearDistance = tempDistance;
							}
						}
					}
				}
			}else if(ModTemp->Type == SPHERE){
				if(ModTemp->S != SphereNull){
					SphereTemp = ModTemp->S;
					if(intersect_sphere(SphereTemp,Q,d,HitPoint) == TRUE){
						tempDistance = sqrt(pow(Q[0] - HitPoint[0],2) + pow(Q[1] - HitPoint[1],2) + pow(Q[2] - HitPoint[2],2));
						if(nearDistance > tempDistance){
							for(i=0;i<3;i++){
								HitFace->color[i] = SphereTemp->color[i];
								N[i] = HitPoint[i] - SphereTemp->center[i];
							}
							normalize(N,N);
							for(i=0;i<3;i++)
								HitFace->normal[i] = N[i];
							HitFace->radius = SphereTemp->radius;
							HitFace->Ks = SphereTemp->Ks;
							HitFace->Kt = SphereTemp->Kt;
							HitFace->reIndex = SphereTemp->reIndex;
							nearDistance = tempDistance;
							result = TRUE;
						}
					}
				}
			}
		}
	}
	return result;
}

double  inner_product(double x[3],double y[3]){
	double temp;
	temp = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
	return temp;
}
int shadow_hit_sphere(struct Sphere *TestSphere,double Q[3],double d[3])
{
	double A,B,C,D;
	double t[2];
	A = 1;
	B = 2*(d[0]*(Q[0] - TestSphere->center[0]) + d[1]*(Q[1] - TestSphere->center[1]) + d[2]*(Q[2] -TestSphere->center[2]));
	C = (pow(Q[0] - TestSphere->center[0],2) + pow(Q[1] - TestSphere->center[1],2) + pow(Q[2] - TestSphere->center[2],2) - pow(TestSphere->radius,2));
	D = pow(B,2) - 4*A*C;

	
	if(D <= 0)
		return FALSE;
	else{
		t[0] = (-B+sqrt(D))/2;
		t[1] = (-B-sqrt(D))/2;
		if(t[0] <=0 || t[1] <=0)
			return FALSE;
		return TRUE;
	}
}
int shadow_hit_poly(struct Poly *TestPoly,double Q[3],double d[3],double HitPoint[3])
{
	double FuncA,FuncB,FuncC,FuncD,t,N_d;
//	float A,B,C,D;
	int i,U,V;
	double vertexUV[3][2],Ruv[2];
	double sTra[3],bTra,temp;

	FuncA = TestPoly->normal[0];
	FuncB = TestPoly->normal[1];
	FuncC = TestPoly->normal[2];
	FuncD = -1*(FuncA * TestPoly->vertex[0][0] + FuncB * TestPoly->vertex[0][1] + FuncC * TestPoly->vertex[0][2]);
	N_d = FuncA*d[0] + FuncB*d[1] + FuncC*d[2];//normal dot d vector
	
	if(N_d >= 0)	// no intersection or inisible
		return FALSE;
	t = -1*(FuncA*Q[0] + FuncB*Q[1] + FuncC*Q[2] + FuncD)/(N_d);
	if(t < 0)
		return FALSE;
	
	if(fabs(TestPoly->normal[0]) > fabs(TestPoly->normal[1]))
		if(fabs(TestPoly->normal[0]) > fabs(TestPoly->normal[2])){
			U = 1;
			V = 2;
		}else{
			U = 0;
			V = 1;
		}
	else
		if(fabs(TestPoly->normal[1]) > fabs(TestPoly->normal[2])){
			U = 0;
			V = 2;
		}else{
			U = 0;
			V = 1;
		}
	
	Ruv[0] = Q[U] + d[U]*t;
	Ruv[1] = Q[V] + d[V]*t;
	for(i=0;i<3;i++){
		vertexUV[i][0] = TestPoly->vertex[i][U];
		vertexUV[i][1] = TestPoly->vertex[i][V];
	}
	bTra = fabs(triangle_area(vertexUV[0],vertexUV[1],vertexUV[2]));
	sTra[0] = fabs(triangle_area(Ruv,vertexUV[1],vertexUV[2]));
	sTra[1] = fabs(triangle_area(vertexUV[0],Ruv,vertexUV[2]));
	sTra[2] = fabs(triangle_area(vertexUV[0],vertexUV[1],Ruv));
	temp =fabs( bTra - sTra[0] - sTra[1] - sTra[2])/bTra;
	
	if(temp <= 0.0001){
		return TRUE;
	}else
		return FALSE;
		

	
}
float shadow(double ray[3],double point[3],struct Object *Obj)
{
	double fatt = 1;
	int i;
	struct Object *ObjTemp;
	struct Model *ModTemp;
	struct Poly *PolyTemp;
	struct Sphere *SphereTemp;
	double HitPoint[3];
	PolyTemp = (struct Poly *)malloc(sizeof(struct Poly));
	ObjTemp = Obj;
	while(Obj->next != ObjectNull){
		Obj = Obj->next;
		ModTemp = Obj->model;
		if(ModTemp != ModelNull){
			if(ModTemp->Type == POLY){
				if(ModTemp->subset != PolyNull)
					for(i = 0; i < ModTemp->count; i++){
						PolyTemp = ModTemp->subset[i];	
						if(shadow_hit_poly(PolyTemp,point,ray,HitPoint) == TRUE){
							fatt = PolyTemp->Kt;
						}
					}
			}else{
				if(ModTemp->S != SphereNull){
					SphereTemp = ModTemp->S;
					if(shadow_hit_sphere(SphereTemp,point,ray) == TRUE){
						fatt = SphereTemp->Kt;
					}
				}
			}
		}
	}

	return fatt;

}

void getR(double N[3],double L[3],double R[3])
{
	int i;
	double _N[3],innerNL,S[3];
	innerNL = inner_product(N,L);
	for(i=0;i<3;i++)
		_N[i] = innerNL*N[i];
	for(i=0;i<3;i++)
		S[i] = _N[i] - L[i];
	for(i=0;i<3;i++)
		R[i] = L[i] + 2*S[i];
	normalize(R,R);
}
void getV(double worldCoord[3],double E[3],double V[3])
{
	V[0] = E[0] - worldCoord[0];
	V[1] = E[1] - worldCoord[1];
	V[2] = E[2] - worldCoord[2];

	normalize(V,V);
}
void get_illumiation_color(unsigned char color[3],double Point[3],double normal[3],double L[3],double lightSet[3])
{
	int i;
	double R[3],V[3],innerNL,innerRV,d;
	
	d=sqrt(pow(Point[0] - lightSet[0],2) + pow(Point[1] - lightSet[1],2) + pow(Point[2] - lightSet[2],2));
	//d:distance of light to surface

	innerNL = inner_product(normal,L);
	getR(normal,L,R);		
	getV(Point,eye,V);//V:the vector from surface to eye
	innerRV = inner_product(R,V);
	if(innerNL < 0)innerNL = 0;
	if(innerNL > 1)innerNL = 1;
	if(innerRV < 0)innerRV = 0;
	if(innerRV > 1)innerRV = 1;
	for(i=0;i<3;i++){
		color[i] = (Kd*diffuse_color[i]*innerNL + 
			Ks*pow(innerRV,shininess)*speculor_color[i])/d;
	}
}
int transmit_ray(double ray[3],double normal[3],double n,double nt,double Tray[3])
{

		int xyz;
	
	for(xyz=0;xyz<3;xyz++){	//realistic ray tracing  P46
	Tray[xyz] = n*(ray[xyz] - normal[xyz]*inner_product(ray,normal))/nt 
				- normal[xyz] * sqrt(1 -  pow(n,2)*(1 - pow(inner_product(ray,normal),2))/pow(nt,2));
	}
	
	normalize(Tray,Tray);
	return TRUE;
/*	int xyz;
	float Ci,SqrtTemp;
	Ci = normal[0]*(-ray[0]) + normal[1]*(-ray[1]) + normal[2]*(-ray[2]);
	SqrtTemp = 1 + reIndex*reIndex*(Ci*Ci - 1);
	if(SqrtTemp < 0)
		return FALSE;
	SqrtTemp = sqrt(SqrtTemp);
	for(xyz=0;xyz<3;xyz++){
		Tray[xyz] = reIndex*ray[xyz] + (reIndex*Ci - SqrtTemp)*normal[xyz];
	}
	return TRUE;
*/
}

unsigned char *RT_trace(double Q[3],double ray[3],double ks,int depth,double now_reIndex)
{
	int rgb,xyz;
	double L[3],fatt;
	unsigned char *color,IluColor[3],*rColor,*tColor;
	struct Poly *HitFace;
	struct Light *lightTemp;
	double HitPoint[3],R[3],nray[3],Tray[3];

	HitFace = (struct Poly *)malloc(sizeof(struct Poly));

	color = (unsigned char *)malloc(sizeof(unsigned char)*3);
	
	color[0] = 0;
	color[1] = 0;
	color[2] = 0;
	if(depth > MAX_LEVEL || ks < MIN_KS )
		return color;	//return z_color
	
	//compute the nearest intersection
	if(get_hit(ray,Q,scene,HitFace,HitPoint) == FALSE)
		return color;	//if no hit, return back ground color
	for(rgb = 0 ; rgb < 3 ;rgb++)
		color[rgb] = HitFace->color[rgb]*Ka;
	/**/
	if(HitFace->Kt>0){
		
		if(transmit_ray(ray,HitFace->normal,now_reIndex,HitFace->reIndex,Tray) == TRUE){
			normalize(Tray,Tray);
			tColor = RT_trace(HitPoint,Tray,ks*HitFace->Kt,depth+1,HitFace->reIndex);
			for(rgb=0;rgb<3;rgb++){	
				if(color[rgb] + tColor[rgb] >255)
					color[rgb] = 255;
				else if(color[rgb] + tColor[rgb] < 0)
					color[rgb] = 0;
				else
					color[rgb] += tColor[rgb];
		
			}
		}
	}
	/*---get light---*/
	lightTemp = lightSourceHead;		
	while(lightTemp->next != lightNull){
		lightTemp = lightTemp->next;
		//check if P is directly lighted
		for(xyz=0;xyz<3;xyz++)
			L[xyz] =lightTemp->lightSet[xyz] - HitPoint[xyz];	
		normalize(L,L);	
		fatt = shadow(L,HitPoint,scene);
		if(fatt >0){	
			get_illumiation_color(IluColor,HitPoint,HitFace->normal,L,lightTemp->lightSet);
			for(rgb=0;rgb<3;rgb++){				
				IluColor[rgb] = IluColor[rgb]*fatt;
				if(color[rgb] + IluColor[rgb] > 255)
					color[rgb] = 255;
				else if(color[rgb] + IluColor[rgb] <0)
					color[rgb] = 0;
				else
					color[rgb] = color[rgb] + IluColor[rgb];
			}		
		}		
	}
	/*get reflection color*/
	if(HitFace->Ks>0){
		for(xyz=0;xyz<3;xyz++)
			nray[xyz] = -ray[xyz];
		getR(HitFace->normal,nray,R);		
		rColor = RT_trace(HitPoint,R,ks*HitFace->Ks,depth+1,now_reIndex);
		for(rgb=0;rgb<3;rgb++){
			if(color[rgb] + rColor[rgb] >255)
				color[rgb] = 255;
			else if(color[rgb] + rColor[rgb] <0)
				color[rgb] = 0;
			else
				color[rgb] = color[rgb] + rColor[rgb];
		
		}
	}
	
	for(rgb=0;rgb<3;rgb++)
		color[rgb] *= ks;
	
	return color;
}
double *prime_ray(double Eye[3],double P[3])
{
	int xyz;
	double *rayTemp;
	rayTemp = (double *)malloc(sizeof(double)*3);
	for(xyz = 0; xyz < 3 ;xyz++)
		rayTemp[xyz] = P[xyz] - Eye[xyz];
	normalize(rayTemp,rayTemp);
	return rayTemp;
}
void ray_tracing(double Eye[3],double width,double height)
{
	int xyz,rgb,y,x;
	double u[3],v[3],w[3],C[3],P[4][3],_t[2][3],PTemp[3][3],dP[3],*ray;
	unsigned char *color;
	for(xyz=0;xyz<3;xyz++){	/*calculate u,v,w vector */
		u[xyz] = ViewMat[0][xyz];
		v[xyz] = ViewMat[1][xyz];
		w[xyz] = ViewMat[2][xyz];
	}

	for(xyz=0;xyz<3;xyz++){
		C[xyz] = Eye[xyz] - n*w[xyz];
		P[0][xyz] = C[xyz] + l*u[xyz] +b*v[xyz];
		P[1][xyz] = C[xyz] + r*u[xyz] +b*v[xyz];
		P[2][xyz] = C[xyz] + r*u[xyz] +t*v[xyz];
		P[3][xyz] = C[xyz] + l*u[xyz] +t*v[xyz];
		_t[0][xyz] = (P[3][xyz]-P[0][xyz])/(height - 1);	/*dT1 = (P3 - P0)/(H - 1)*/
		_t[1][xyz] = (P[2][xyz]-P[1][xyz])/(height - 1);	/*dT2 = (P2 - P1)/(H - 1)*/
		
		PTemp[1][xyz] = P[0][xyz];	/*P' = P0*/
		PTemp[2][xyz] = P[1][xyz];	/*P'' = P1*/	
	}

	for(y=0;y<height;y++){		/*calculate each rays*/
		for(xyz=0;xyz<3;xyz++){
			dP[xyz] = (PTemp[2][xyz] - PTemp[1][xyz])/(width - 1);	/*dP = (P'' - P')/(W-1)*/
			PTemp[0][xyz] = PTemp[1][xyz];						/*P = P'*/
		}
		for(x=0;x<width;x++){
			ray = prime_ray(Eye,PTemp[0]);		/*get ray direction*/
			color = RT_trace(Eye,ray,1,1,1);	/*calculate ray*/
			for(rgb=0;rgb<3;rgb++)
				color_buffer[y][x][rgb] = color[rgb];
			free(ray);
			free(color);
			for(xyz = 0;xyz<3;xyz++)
				PTemp[0][xyz] += dP[xyz];
		}
		for(xyz = 0;xyz<3;xyz++){
			PTemp[1][xyz] += _t[0][xyz];
			PTemp[2][xyz] += _t[1][xyz];
		}
	}
}

static void draw_image(unsigned char img[][WSIZE][3],int xsize,int ysize)
{
	int i,j;
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(xsize/2,ysize/2,100,xsize/2,ysize/2,-1.0,0,1,0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-xsize/2,xsize/2,-ysize/2,ysize/2,-100.0,100.0);


	glBegin(GL_POINTS);
	for(i=0;i<ysize;i++)
		for(j=0;j<xsize;j++){
			glColor3ub(img[i][j][0],img[i][j][1],img[i][j][2]);		
			glVertex2f((double)j,(double)i);
		}
	

	glEnd();
	
}
void my_display_func(void)
{
	init_func();
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	LookAt(eye[0],eye[1],eye[2],0,0,0,0,1,0);	
	projection(-1,-1,1,1,1,100);

	PushM();
	Rotate(EyeRota,0,1,0);
	draw();
	PopM();	
	ray_tracing(eye,WSIZE,WSIZE);
	draw_image(color_buffer,WSIZE,WSIZE);

	glutSwapBuffers();
	glFlush();
}
void my_reshape_func(int w,int h)
{
	
	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    /* setup the projection matrix, such that the rectangle is correctly
       drawn
       */
    if (w <= h) 
      glOrtho (-50.0, 50.0, -50.0*(GLfloat)h/(GLfloat)w, 
	       50.0*(GLfloat)h/(GLfloat)w, -1.0, 1.0);
    else 
      glOrtho (-50.0*(GLfloat)w/(GLfloat)h, 
	       50.0*(GLfloat)w/(GLfloat)h, -50.0, 50.0, -1.0, 1.0);

    glViewport(0, 0, w, h);

    /*--- Set current mtx = modelview matrix ---*/
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity ();       /*--- clear modelview matrix ---*/

}
void Exit()
{
	int PolyCount;
	struct Poly *Ptemp;
	struct Sphere *Stemp;
	struct Model *Mtemp;
	struct Object *Otemp,*Odelete;
	Otemp = scene;
	while(Otemp!=ObjectNull){
		Mtemp = Otemp->model;
		if(Mtemp!=ModelNull){
			if(Mtemp->Type == POLY){
				for(PolyCount = 0;PolyCount<Mtemp->count;PolyCount++){
					Ptemp = Mtemp->subset[PolyCount];
					free(Ptemp);
				}
			}else{
				Stemp = Mtemp->S;
				free(Stemp);
			}
			free(Mtemp);
		}
		Odelete = Otemp;
		Otemp = Otemp->next;
		free(Odelete);
	}
	
	free(PolyNull);
	free(ModelNull);
	free(ObjectNull);
	exit(0);
}
void my_keyboard_func(unsigned char key,int ix,int iy)
{
	if(key == 'q' || key == 'Q')
		Exit();
	if(key == '0'){
		EyeRota += 15;	
		my_display_func();
	}	
}
void my_mouse_func(int button,int state,int x,int y)
{


}
void my_motion_func(int x,int y)
{
}
void my_idle_func(void)
{
}

void main(int argc,char **argv)
{
	glutInit(&argc,argv);	 /*---Make a session with windows system---*/

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); /*﹚竡display篈*/
	glutInitWindowSize(WSIZE,WSIZE);/*window size (width,height)*/
	glutCreateWindow("Ray Tracing");

	init_func();

	glutDisplayFunc(my_display_func); /* display event callback func */
	glutReshapeFunc(my_reshape_func); /* reshape event callback func */
	glutKeyboardFunc(my_keyboard_func);/* keyboard event callback func */
	glutMouseFunc(my_mouse_func);  /* Mouse Button ecent callback func */
	glutMotionFunc(my_motion_func);/* Mouse motion event callback func */
	glutIdleFunc(my_idle_func);

	glutMainLoop(); 
}