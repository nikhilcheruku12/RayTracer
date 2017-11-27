/*
CSCI 480
Assignment 3 Raytracer

Name: <Nikhil Cherukuri>
*/

#include <stdlib.h>
#include <stdio.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include <math.h> 
#include <iostream>
#include <fstream>
#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//#define WIDTH 320
//#define HEIGHT 240

//the field of view of the camera
#define fov 60.0
#define SCREEN_DISTANCE 1

#define PI 3.14159265

unsigned char buffer[HEIGHT][WIDTH][3];
double SCREEN_WIDTH;
double SCREEN_HEIGHT;

std::ofstream outfile;
double SCREEN_HEIGHT_INCR;
double SCREEN_WIDTH_INCR;


struct point {
   double x;
   double y;
   double z;
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Ray
{
  double* origin;
  double* direction;
};

Ray rays[HEIGHT][WIDTH];
typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene()
{
  /*glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
  glLoadIdentity();             // Reset
   // Enable perspective projection with fovy, aspect, zNear and zFar
  gluPerspective(60.0f, (GLfloat)WIDTH/(GLfloat)HEIGHT, 0.01f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  gluLookAt(0,0,0,
            0,0,-1,
            0,1,0);*/
  unsigned int x,y;

  
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      //plot_pixel(x,y,x%256,y%256,(x+y)%256);
      plot_pixel(x,y,buffer[y][x][0],buffer[y][x][1],buffer[y][x][2]);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

point make_unit (point p){
    float divisor = sqrt (pow(p.x,2) + pow(p.y,2) + pow(p.z,2));
    if(divisor == 0)
    {
      p.x = 0; p.y = 0; p.z = 0;
      return p;
    }
    p.x = p.x / divisor;
    p.y = p.y / divisor;
    p.z = p.z / divisor;
    return p;
  }

double distance (point p1, point p2){
  return sqrt( pow(p1.x - p2.x,2) + pow(p1.y - p2.y,2) + pow(p1.z - p2.z,2));
}

point getPoint (double* d){
  point p; p.x = d[0]; p.y = d[1]; p.z = d[2];
}
double* makeUnit (double* p){
    double divisor = sqrt (pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
    if(divisor == 0)
    {
      p[0] = 0; p[1] = 0; p[2] = 0;
      return p;
    }
    double p1[3];
    std::cout << "p before = " << "{ " << p[0] << ", " << p[1] << ", " << p[2] << "}" << std::endl;
    p1[0] = p[0] / divisor;
    p1[1] = p[1] / divisor;
    p1[2] = p[2] / divisor;
    std::cout << "p after = " << "{ " << p1[0] << ", " << p1[1] << ", " << p1[2] << "}" << std::endl;
    //return p1;
    return p1;
  }

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void printRay (Ray ray){
   //outfile<< "origin = " << "{" << ray.origin[0] << ", " << ray.origin[1] << ", " << ray.origin[2] << " }" << std::endl;
  // outfile << "dir = " << "{" << ray.direction[0] << ", " << ray.direction[1] << ", " << ray.direction[2] << " }" << std::endl;
}

void printRay2 (Ray ray){
   //outfile<< "origin = " << "{" << ray.origin[0] << ", " << ray.origin[1] << ", " << ray.origin[2] << " }" << std::endl;
   //std::cout << "dir = " << "{" << ray.direction[0] << ", " << ray.direction[1] << ", " << ray.direction[2] << " }" << std::endl;
}

void crossProduct (double* u, double* v, double* product){
     //double product[3];
    // u_2v_3-u_3v_2,\:u_3v_1-u_1v_3,\:u_1v_2-u_2v_1
     //u1*v2 -  u2*v1 ; u2*v0 - u0*v2; u0*v1 - u1*v0
      product[0] = u[1] * v[2] - u[2] * v[1];
      product[1] = u[2] * v[0] - u[0] * v[2];
      product[2] = u[0] * v[1] - u[1] * v[0];
      
  }

double innerProduct (double* u, double* v){
     double product;
      product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
      return product;
  }

point crossProduct (point u, point v){
     point product;
      product.x = u.y * v.z - v.y * u.z;
      product.y = v.x * u.z - u.x * v.z;
      product.z = u.x * v.y - v.x * u.y;
      return make_unit(product);
  }

double dotProduct (point u, point v){
     double product;
      product = u.x*v.x + u.y*v.y + u.z*v.z;
      return product;
  }

point multiplyPoint (point p1, double t){
  point p; 
  p.x = p1.x * t;
  p.y = p1.y * t;
  p.z = p1.z * t;
  return p;
}

point addPoint (point p1, point p2){
  point p; 
  p.x = p1.x + p2.x;
  p.y = p1.y + p2.y;
  p.z = p1.z + p2.z;
  return p;
}


void vector (double* a, double* b, double* c){
  a[0] = b[0] - c[0]; a[1] = b[1] - c[1]; a[2] = b[2] - c[2];
  return;
}

point getEdge  (point p1, point p2){
  point p; 
  p.x = p1.x - p2.x;
  p.y = p1.y - p2.y;
  p.z = p1.z - p2.z;
  return p;
}

/*bool RayIntersectsTriangle(point rayOrigin, 
                           point rayVector, 
                           Triangle inTriangle,
                           point& outIntersectionPoint)
{
    const float EPSILON = 0.00001; 
    point vertex0 = getPoint(inTriangle.v[0].position);
    point vertex1 = getPoint(inTriangle.v[1].position);  
    point vertex2 = getPoint(inTriangle.v[2].position);
    point edge1, edge2, h, s, q;
    double a,f,u,v;
    edge1 = getEdge(vertex1 , vertex0);
    edge2 = getEdge(vertex2 , vertex0);
    h = crossProduct(rayVector,edge2);
    a = dotProduct(edge1,h);
    if (a > -EPSILON && a < EPSILON)
        return false;
    f = 1/a;
    s =  getEdge(rayOrigin, vertex0);
    u = f * (dotProduct(s,h));
    if (u < 0.0 || u > 1.0)
        return false;
    q = crossProduct(s,edge1);
    v =  f * dotProduct(rayVector,q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double  t =  f * dotProduct(edge2,q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = addPoint(rayOrigin, multiplyPoint(rayVector,t));
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}*/

void printDoubles (double* d, std::string Name){
    outfile << Name << " = {" << d[0] << " , " << d[1] << " , " << d[2] << " }" << std::endl;
}
bool trinagle_intr (Triangle triangle, Ray ray, bool print){
  double* v0 = triangle.v[0].position;
  double* v1 = triangle.v[1].position;
  double* v2 = triangle.v[2].position;
  if(print){
    printDoubles(v0, "v0");
    printDoubles(v1,"v1");
    printDoubles(v2,"v2");
  }
  
  double e1[3],e2[3];
  vector(e1,v1,v0);
  vector(e2,v2,v0);

  if(print){
    printDoubles(e1,"e1");
    printDoubles(e2,"e2");
  }
  

  //double* N = crossProduct(e1,e2); 
  double N[3];
  crossProduct(e1,e2,N);
  //double* N = crossProduct(e2,e1); 
  //double D = innerProduct(N,v0) * -1; 
  //double D = innerProduct(N,v0) ; 
  if(print)
   {
    printDoubles(N, "N");
   } 

   double* temp = makeUnit(N);
   N[0] = temp[0]; N[1] = temp[1]; N[2] = temp[2];
   if(print){
    printDoubles(N, "N unit");
   }

   double D = innerProduct(N,v0) ; 
  //float t = - (dot(N, orig) + D) / dot(N, dir); 
  if( innerProduct(N, ray.direction) == 0)
    return false;
  if (print){
    outfile << "D = " << D << std::endl;
    printDoubles(ray.direction, "direction");
    outfile <<  " innerProduct(N, ray.direction) = " << innerProduct(N, ray.direction) << std::endl;
  }
  double t =  (D) / ((double)innerProduct(N, ray.direction));
  if(print)
  outfile << "t = " << t << std::endl;
  if(t < 0)
    return false;
  


  /*point intr_point = multiplyPoint(getPoint(ray.direction),t);
  double p[3]; p[0] = intr_point.x; p[1] = intr_point.y; p[2] = intr_point.z;*/
  double p[3]; p[0] = t*ray.direction[0]; p[1] = t*ray.direction[1]; p[2] = t*ray.direction[2];
  if(print)
  printDoubles(p, "intr_point");
  double c0[3] ,c1[3] ,c2[3];
  double edg0[3], edg1[3], edg2[3];
  vector(c0,p,v0); vector(c1,p,v1); vector(c2,p,v2);
  vector(edg0,v1,v0); vector(edg1,v2,v1); vector(edg2,v0,v2);

  double temp0[3], temp1[3], temp2[3];
  crossProduct(edg0,c0, temp0);
  crossProduct(edg1,c1, temp1);
  crossProduct(edg2,c2, temp2);
  if (innerProduct(N, temp0) > 0 && 
    innerProduct(N, temp1) > 0 && 
    innerProduct(N, temp2) > 0) return true;

    return false;
}
/*bool triangle_intersection( Triangle triangle, Ray ray) {

  double* p = ray.origin;
  double* d = ray.direction;

  double* v0 = triangle.v[0].position;
  double* v1 = triangle.v[1].position;
  double* v2 = triangle.v[2].position;
  double e1[3],e2[3],s[3];
  double a,f,u,v;
  vector(e1,v1,v0);
  vector(e2,v2,v0);

  double* h = crossProduct(d,e2);
   a = innerProduct(e1,h);

  if (a > -0.00001 && a < 0.00001)
    return(false);

  f = 1/a;
  vector(s,p,v0);
  u = f * (innerProduct(s,h));

  if (u < 0.0 || u > 1.0)
    return(false);

  double* q = crossProduct(s,e1);
  v = f * innerProduct(d,q);

  if (v < 0.0 || u + v > 1.0)
    return(false);

  // at this stage we can compute t to find out where
  // the intersection point is on the line
  double t = f * innerProduct(e2,q);

  if (t > 0.00001) // ray intersection
    return(true);

  else // this means that there is a line intersection
     // but not a ray intersection
     return (false);

}*/


bool sphere_intersection (Sphere sphere, Ray ray){

  double xo = ray.origin[0]; double yo = ray.origin[1]; double zo = ray.origin[2];
  double xd = ray.direction[0]; double yd = ray.direction[1]; double zd = ray.direction[2];
  double xc = sphere.position[0]; double yc = sphere.position[1]; double zc = sphere.position[2];

  double r = sphere.radius;

  double a = (pow(xd,2) + pow(yd,2) + pow(zd,2));
  double b = 2*(xd*(xo-xc) + yd*(yo-yc) + zd*(zo - zc));
  double c = pow ((xo-xc),2) + pow ((yo-yc),2) + pow ((zo-zc),2) - pow(r,2);

  double d = pow(b,2) - 4*(a)*(c);
 // outfile << "d = " << d << " xd = " << xd << " yd = " << yd << " zd = " << zd << std::endl;
  if(d < 0)
    return false;
  else
    return true;
}

/*point sphere_intersection (Sphere sphere, Ray ray){

  double xo = ray.origin[0]; double yo = ray.origin[1]; double zo = ray.origin[2];
  double xd = ray.direction[0]; double yd = ray.direction[1]; double zd = ray.direction[2];
  double xc = sphere.position[0]; double yc = sphere.position[1]; double zc = sphere.position[2];

  double r = sphere.radius;

  double a = (pow(xd,2) + pow(yd,2) + pow(zd,2));
  double b = 2*(xd*(xo-xc) + yd*(yo-yc) + zd*(zo - zc));
  double c = pow ((xo-xc),2) + pow ((yo-yc),2) + pow ((zo-zc),2) - pow(r,2);

  double d = pow(b,2) - 4*(a)*(c);
  outfile << "d = " << d << " xd = " << xd << " yd = " << yd << " zd = " << zd << std::endl;
  if(d < 0)
    return NULL;
  else if (d == 0)
  {
    double t = -b/2*(a);
    point p; p.x = xo + xd*t; p.y = yo + yd*t; p.z = zo + zd*t;
    return p;
  }

  else{
    double t1 = (-b+d)/(2*(a));
    double t2 = (-b-d)/(2*(a));
    point p;
    point p1; p.x = xo + xd*t1; p.y = yo + yd*t1; p.z = zo + zd*t1;
    point p2; p.x = xo + xd*t2; p.y = yo + yd*t2; p.z = zo + zd*t2;

    point rayDir; rayDir.x = ray.direction[0]; rayDir.y = ray.direction[1]; rayDir.z = ray.direction[2];
    double distance1 = distance(p1,rayDir);
    double distance2 = distance(p2,,rayDir);

    if(distance1 < distance2){
      p = p1;
    } else{
      p = p2;
    }
    return p;
  }
}*/

void raySetup (){
  //double angle = fov/2;
  double angle = fov/2;
  double tangent  =  tan ( angle * PI / 180.0 ); 
  SCREEN_WIDTH = tangent * SCREEN_DISTANCE * 2;
  SCREEN_HEIGHT = SCREEN_WIDTH * ((double)((double)HEIGHT/(double)WIDTH));
  //SCREEN_HEIGHT = SCREEN_WIDTH;
  SCREEN_WIDTH_INCR = SCREEN_WIDTH / (double) WIDTH;
  SCREEN_HEIGHT_INCR = SCREEN_HEIGHT /(double) HEIGHT;

  
  double left_most = SCREEN_WIDTH/2 * -1;
  double bottom_most = SCREEN_HEIGHT/2 * -1;

  for (int i = 0; i < HEIGHT; i++ ){
    for (int j = 0; j < WIDTH; j++){
      double x_cord = j*SCREEN_WIDTH_INCR + left_most;
      double y_cord = i*SCREEN_HEIGHT_INCR + bottom_most;

      //double x_cord = SCREEN_WIDTH/2 - (j*SCREEN_WIDTH_INCR);
      //double y_cord = SCREEN_HEIGHT/2 - (i*SCREEN_HEIGHT_INCR);

      Ray ray;
      //double temp[3] = {0,0,0};
      ray.origin = new double[3]; ray.origin[0] = 0; ray.origin[1] =0; ray.origin[2] = 0;
      //double temp2[3] = {x_cord, y_cord, SCREEN_DISTANCE*-1};
      //ray.direction = temp2;
      ray.direction = new double[3]; ray.direction[0] = x_cord; ray.direction[1] = y_cord; ray.direction[2] = SCREEN_DISTANCE * -1;
      //ray.direction = makeUnit(ray.direction);
      double* temp = makeUnit(ray.direction);
      ray.direction[0] = temp[0]; ray.direction[1] = temp[1]; ray.direction[2] = temp[2];
      rays[i][j] = ray;
      bool temp3 = sphere_intersection(spheres[0],rays[i][j]);
      printRay2(rays[i][j]);
       //outfile << "intersects = " << temp3 << std::endl;

    }
  }

  for (int i = 0; i < HEIGHT; i++ ){
    for (int j = 0; j < WIDTH; j++){
      //outfile << " i = " << i << " j = " << j << " ";
      printRay(rays[i][j]);
       //outfile << "intersects = " << temp3 << std::endl;
    }
  }

  std::cout << "angle = " << angle << std::endl;
  std::cout << "tangent = " << tangent << std::endl;
  std::cout << "SCREEN_WIDTH = " << SCREEN_WIDTH << std::endl;
  std::cout << "SCREEN_HEIGHT = " << SCREEN_HEIGHT << std::endl;
  std::cout << "SCREEN_WIDTH_INCR = " << SCREEN_WIDTH_INCR << std::endl;
  std::cout << "SCREEN_HEIGHT_INCR = " << SCREEN_HEIGHT_INCR << std::endl;

  for (int i = 0; i < num_spheres; i++){
    //if(sphere_intersection)

    for (int j = 0; j < HEIGHT; j++){
      for (int k = 0; k < WIDTH; k++){
         bool intersection = sphere_intersection(spheres[i],rays[j][k]);
         if(intersection){
           buffer[j][k][0] = 255; buffer[j][k][1] = 255; buffer[j][k][2] = 255;  
         } else{
           //buffer[j][k][0] = 0; buffer[j][k][1] = 0; buffer[j][k][2] = 0;
         }
         
      }
    }
  }

  for (int i = 0; i < num_triangles; i++){
    //if(sphere_intersection)
    outfile << "Triangle number = " << i << std::endl;
    for (int j = 0; j < HEIGHT; j++){
      for (int k = 0; k < WIDTH; k++){
         //bool intersection = triangle_intersection(triangles[i],rays[j][k]);
         // point p;
         /*bool intersection = RayIntersectsTriangle(getPoint(rays[j][k].origin),
                                                   getPoint(rays[j][k].direction),
                                                   triangles[i], p );*/
         //outfile << "k = " << k << " j = " << j << std::endl;
          bool intersection;
          if(k >= 67 && k <= 272 && j >= 31 && j <= 235){
            outfile << "k = " << k << " j = " << j << std::endl;
            intersection = trinagle_intr(triangles[i], rays[j][k], true);
          }
           

          else
            intersection = trinagle_intr(triangles[i], rays[j][k], false);
         if(intersection){
           buffer[j][k][0] = 255; buffer[j][k][1] = 255; buffer[j][k][2] = 255;  
         } else{
           //buffer[j][k][0] = 0; buffer[j][k][1] = 0; buffer[j][k][2] = 0;
         }
         
      }
    }
  }
  //printf(" SCREEN_WIDTH = %s , SCREEN_HEIGHT = %s\n", SCREEN_WIDTH, SCREEN_HEIGHT); 
  outfile.close();
}

double* getDir(int ix, int iy){
  double direction[3];
  direction[2] = SCREEN_DISTANCE*-1;
  return direction;
}
void renderPixel(int ix, int iy)
{
    Ray ray;
    bool objFound = false;  // Not used. Only used for reflect traces.
    double temp[3];
    temp[0] = 0;  temp[1] = 0; temp[2] = 0;
    ray.origin = temp;
    ray.direction = getDir(ix, iy);

    //vec4 color = trace(ray, objFound);
    //setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < HEIGHT; iy++)
        for (int ix = 0; ix < WIDTH; ix++)
            renderPixel(ix, iy);
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  /*glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
  glLoadIdentity();             // Reset
   // Enable perspective projection with fovy, aspect, zNear and zFar
  gluPerspective(60.0f, aspect, 0.01f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);*/

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
  outfile.open(("PersonalData.txt"));
  raySetup();

}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
