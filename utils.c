/*
   utils.c - F.J. Estrada, Dec. 9, 2010

   Utilities for the ray tracer. You will need to complete
   some of the functions in this file. Look for the sections
   marked "TO DO". Be sure to read the rest of the file and
   understand how the entire code works.
*/

#include "utils.h"

// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4]={{1.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0}};

/////////////////////////////////////////////
// Primitive data structure section
/////////////////////////////////////////////
struct point3D *newPoint(double px, double py, double pz)
{
 // Allocate a new point structure, initialize it to
 // the specified coordinates, and return a pointer
 // to it.

 struct point3D *pt=(struct point3D *)calloc(1,sizeof(struct point3D));
 if (!pt) fprintf(stderr,"Out of memory allocating point structure!\n");
 else
 {
  pt->px=px;
  pt->py=py;
  pt->pz=pz;
  pt->pw=1.0;
 }
 return(pt);
}

struct pointLS *newPLS(struct point3D *p0, double r, double g, double b)
{
 // Allocate a new point light sourse structure. Initialize the light
 // source to the specified RGB colour
 // Note that this is a point light source in that it is a single point
 // in space, if you also want a uniform direction for light over the
 // scene (a so-called directional light) you need to place the
 // light source really far away.

 struct pointLS *ls=(struct pointLS *)calloc(1,sizeof(struct pointLS));
 if (!ls) fprintf(stderr,"Out of memory allocating light source!\n");
 else
 {
  memcpy(&ls->p0,p0,sizeof(struct point3D));	// Copy light source location

  ls->col.R=r;					// Store light source colour and
  ls->col.G=g;					// intensity
  ls->col.B=b;
 }
 return(ls);
}

/////////////////////////////////////////////
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj)
{
 // Transforms a ray using the inverse transform for the specified object. This is so that we can
 // use the intersection test for the canonical object. Note that this has to be done carefully!

 ///////////////////////////////////////////
 // TO DO: Complete this function
 ///////////////////////////////////////////
  memcpy(&ray_transformed->p0, &ray_orig->p0, sizeof(struct point3D));
  memcpy(&ray_transformed->d, &ray_orig->d, sizeof(struct point3D));
  
  matVecMult(obj->Tinv, &ray_transformed->p0);
  matVecMult(obj->Tinv, &ray_transformed->d);

  ray_transformed->rayPos=&rayPosition;
}

inline void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj)
{
 // Computes the normal at an affinely transformed point given the original normal and the
 // object's inverse transformation. From the notes:
 // n_transformed=A^-T*n normalized.

 ///////////////////////////////////////////
 // TO DO: ComTe this function
 ///////////////////////////////////////////
 double C[4][4];
 int i,j,k;

 memset(C,0,16*sizeof(double));
 for (i=0;i<4;i++)
  for (j=0;j<4;j++)
     C[i][j]=obj->Tinv[j][i];

 memcpy(n_transformed,n_orig,sizeof(struct point3D));
 normalize(n_transformed);

 matVecMult(C, n_transformed);

}

/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////
struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
 // Intialize a new plane with the specified parameters:
 // ra, rd, rs, rg - Albedos for the components of the Phong model
 // r, g, b, - Colour for this plane
 // alpha - Transparency, must be set to 1 unless you are doing refraction
 // r_index - Refraction index if you are doing refraction.
 // shiny - Exponent for the specular component of the Phong model
 //
 // The plane is defined by the following vertices (CCW)
 // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
 // With normal vector (0,0,1) (i.e. parallel to the XY plane)

 struct object3D *plane=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!plane) fprintf(stderr,"Unable to allocate new plane, out of memory!\n");
 else
 {
  plane->alb.ra=ra;
  plane->alb.rd=rd;
  plane->alb.rs=rs;
  plane->alb.rg=rg;
  plane->col.R=r;
  plane->col.G=g;
  plane->col.B=b;
  plane->alpha=alpha;
  plane->r_index=r_index;
  plane->shinyness=shiny;
  plane->intersect=&planeIntersect;
  plane->texImg=NULL;
  memcpy(&plane->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&plane->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  plane->textureMap=&texMap;
  plane->frontAndBack=1;
  plane->isLightSource=0;
 }
 return(plane);
}

struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
 // Intialize a new sphere with the specified parameters:
 // ra, rd, rs, rg - Albedos for the components of the Phong model
 // r, g, b, - Colour for this plane
 // alpha - Transparency, must be set to 1 unless you are doing refraction
 // r_index - Refraction index if you are doing refraction.
 // shiny -Exponent for the specular component of the Phong model
 //
 // This is assumed to represent a unit sphere centered at the origin.
 //

 struct object3D *sphere=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!sphere) fprintf(stderr,"Unable to allocate new sphere, out of memory!\n");
 else
 {
  sphere->alb.ra=ra;
  sphere->alb.rd=rd;
  sphere->alb.rs=rs;
  sphere->alb.rg=rg;
  sphere->col.R=r;
  sphere->col.G=g;
  sphere->col.B=b;
  sphere->alpha=alpha;
  sphere->r_index=r_index;
  sphere->shinyness=shiny;
  sphere->intersect=&sphereIntersect;
  sphere->texImg=NULL;
  memcpy(&sphere->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&sphere->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  sphere->textureMap=&texMap;
  sphere->frontAndBack=0;
  sphere->isLightSource=0;
 }
 return(sphere);
}

struct object3D *newCylinder(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
 // Intialize a new cylinder with the specified parameters:
 // ra, rd, rs, rg - Albedos for the components of the Phong model
 // r, g, b, - Colour for this plane
 // alpha - Transparency, must be set to 1 unless you are doing refraction
 // r_index - Refraction index if you are doing refraction.
 // shiny -Exponent for the specular component of the Phong model
 //
 // This is assumed to represent a unit sphere centered at the origin.
 //

 struct object3D *cylinder=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!cylinder) fprintf(stderr,"Unable to allocate new sphere, out of memory!\n");
 else
 {
  cylinder->alb.ra=ra;
  cylinder->alb.rd=rd;
  cylinder->alb.rs=rs;
  cylinder->alb.rg=rg;
  cylinder->col.R=r;
  cylinder->col.G=g;
  cylinder->col.B=b;
  cylinder->alpha=alpha;
  cylinder->r_index=r_index;
  cylinder->shinyness=shiny;
  cylinder->intersect=&cylinderIntersect;
  cylinder->texImg=NULL;
  memcpy(&cylinder->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&cylinder->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  cylinder->textureMap=&texMap;
  cylinder->frontAndBack=0;
  cylinder->isLightSource=0;
 }
 return(cylinder);
}



///////////////////////////////////////////////////////////////////////////////////////
// TO DO:
//	Complete the functions that compute intersections for the canonical plane
//      and canonical sphere with a given ray. This is the most fundamental component
//      of the raytracer.
///////////////////////////////////////////////////////////////////////////////////////
void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical plane.

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
 struct ray3D * origin;
 struct point3D *tn, *new_ray_d, *new_ray_p;
 new_ray_d = newPoint(0 ,0 ,0);
 new_ray_p = newPoint(0, 0, 0);
 new_ray_d->pw = 0;
 origin = newRay(new_ray_p, new_ray_d);
 rayTransform(ray, origin, plane);
 *lambda = -1;
 if(origin->d.pz == 0 || plane->isLightSource ==1){
  free(new_ray_d);
 free(new_ray_p);
 free(origin);
  return;
 }
 double t = (-origin->p0.pz) / origin->d.pz;
 if(t < 0 || origin->d.pz == 0){
 free(new_ray_d);
 free(new_ray_p);
 free(origin);
  return;
 }

 rayPosition(origin, t, p);
 tn = newPoint(0, 0, 1); 

 if (plane->texImg != NULL){
      *a = max(0, -1*max(-1, -(p->px+1)/2));
      *b = max(0, -1*max(-1, -(p->py+1)/2));
    }

 if(p->px >= -1 && p->px <= 1 && p->py >= -1 && p->py <= 1){

  *lambda = t;
  normalTransform(tn, n, plane); 
  matVecMult(plane->T, p);
 }
 free(new_ray_d);
 free(new_ray_p);
 free(origin);
 free(tn);

}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical sphere.

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
  struct ray3D * origin;
  struct point3D *tn, *new_ray_d, *new_ray_p;
 new_ray_d = newPoint(0 ,0 ,0);
 new_ray_p = newPoint(0, 0, 0);
 new_ray_p->pw = 1;
 new_ray_d->pw = 0;
 origin = newRay(new_ray_p, new_ray_d);
 rayTransform(ray, origin, sphere);
 *lambda = -1;

 double t, t1, t2, A, B, C;

 A=dot(&origin->d, &origin->d);
 B = 2* dot(&origin->d, &origin->p0);
 C = dot(&origin->p0, &origin->p0) -1;
 if(B*B - 4 *A*C < 0){
  free(new_ray_d);
 free(new_ray_p);
 free(origin);
  return;
 }
 t1 = (-B - sqrt(B*B - 4 *A*C))/(2*A);
 t2 = (-B + sqrt(B*B - 4 *A*C))/(2*A);

 if(t1 < t2 && t1 > 0){
   t = t1;
 }
 else{
   t = t2;
 }

 if (t <= 0){
  free(new_ray_d);
  free(new_ray_p);
  free(origin);
  return;
 }
 rayPosition(origin, t, p);
 
 tn = newPoint(-p->px, -p->py, -p->pz);

 *lambda = t;
 normalTransform(tn, n, sphere);
 matVecMult(sphere->T, p);

 free(new_ray_d);
 free(new_ray_p);
 free(origin);
 free(tn);
}

void cylinderIntersect(struct object3D *cylinder, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical sphere.

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
 struct ray3D * origin;
 struct point3D *tn, *new_ray_d, *new_ray_p;
 new_ray_d = newPoint(0 ,0 ,0);
 new_ray_p = newPoint(0, 0, 0);
 new_ray_p->pw = 1;
 new_ray_d->pw = 0;
 origin = newRay(new_ray_p, new_ray_d);
 rayTransform(ray, origin, cylinder);
 *lambda = -1;

 double t, t1, t2, A, B, C, z1, z2, z3;
 double p0_z = origin->p0.pz;
 double d_z = origin->d.pz;
 origin->p0.pz = 0;
 origin->d.pz = 0;
 A=dot(&origin->d, &origin->d);
 B = 2* dot(&origin->d, &origin->p0);
 C = dot(&origin->p0, &origin->p0) -1;
 if(B*B - 4 *A*C < 0){
  free(new_ray_d);
 free(new_ray_p);
 free(origin);
  return;
 }
 t1 = (-B - sqrt(B*B - 4 *A*C))/(2*A);
 t2 = (-B + sqrt(B*B - 4 *A*C))/(2*A);


 z1 = p0_z +(t1*d_z);
 z2 = p0_z +(t2*d_z);
 t =-1;

 if(z1 <= 1 && z1 >= 0){
  t = z1;	
 }
 if(z2 <= 1 && z2 >= 0 && z2 > t){
  t = z2;
 }  

 if (t <= 0){
  free(new_ray_d);
  free(new_ray_p);
  free(origin);
  return;
 }

 rayPosition(origin, t, p);

 tn = newPoint(-p->px, -p->py, -p->pz);

 *lambda = t;
 normalTransform(tn, n, cylinder);
 matVecMult(cylinder->T, p);
 
 free(new_ray_d);
 free(new_ray_p);
 free(origin);
 free(tn);
}


void loadTexture(struct object3D *o, const char *filename)
{
 // Load a texture image from file and assign it to the
 // specified object
 if (o!=NULL)
 {
  if (o->texImg!=NULL)	// We have previously loaded a texture
  {			// for this object, need to de-allocate it
   if (o->texImg->rgbdata!=NULL) free(o->texImg->rgbdata);
   free(o->texImg);
  }
  o->texImg=readPPMimage(filename);	// Allocate new texture

 }
}

void texMap(struct image *img, double a, double b, double *R, double *G, double *B)
{
 /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
   a given object.

  The colour is returned in R, G, B. Uses bi-linear interpolation
  to determine texture colour.
 */

 //////////////////////////////////////////////////
 // TO DO (Assignment 4 only):
 //
 //  Complete this function to return the colour
 // of the texture image at the specified texture
 // coordinates. Your code should use bi-linear
 // interpolation to obtain the texture colour.
 //////////////////////////////////////////////////
  int i, j;
  double up, vp;  
  double *rgbIm;

  rgbIm=(double *)img->rgbdata;

  i = floor(img->sx*a);
  j = floor(img->sy*b);

  up = img->sx*a - i;
  vp = img->sy*b - j;

  *(R) = (1-up)*(1-vp)*rgbIm[(j*img->sx+i)*3] + up*(1-vp)*rgbIm[((j+1)*img->sx+i)*3] + (1-up)*vp*rgbIm[(j*img->sx+i+1)*3] + up*vp*rgbIm[((j+1)*img->sx+i+1)*3];
  
  *(G) = (1-up)*(1-vp)*rgbIm[(j*img->sx+i)*3 + 1] + up*(1-vp)*rgbIm[((j+1)*img->sx+i)*3 +1] + (1-up)*vp*rgbIm[(j*img->sx+i+1)*3 + 1] + up*vp*rgbIm[((j+1)*img->sx+i+1)*3 +1];
  
  *(B) = (1-up)*(1-vp)*rgbIm[(j*img->sx+i)*3 + 2] + up*(1-vp)*rgbIm[((j+1)*img->sx+i)*3 + 2] +  (1-up)*vp*rgbIm[(j*img->sx+i+1)*3 + 2] + up*vp*rgbIm[((j+1)*img->sx+i+1)*3 + 2];

 return;
}

void insertObject(struct object3D *o, struct object3D **list)
{
 if (o==NULL) return;
 // Inserts an object into the object list.
 if (*(list)==NULL)
 {
  *(list)=o;
  (*(list))->next=NULL;
 }
 else
 {
  o->next=(*(list))->next;
  (*(list))->next=o;
 }
}

void insertPLS(struct pointLS *l, struct pointLS **list)
{
 if (l==NULL) return;
 // Inserts a light source into the list of light sources
 if (*(list)==NULL)
 {
  *(list)=l;
  (*(list))->next=NULL;
 }
 else
 {
  l->next=(*(list))->next;
  (*(list))->next=l;
 }

}

void addAreaLight(float sx, float sy, float nx, float ny, float nz,\
                  float tx, float ty, float tz, int lx, int ly,\
                  float r, float g, float b, struct object3D **o_list, struct pointLS **l_list)
{
 /*
   This function sets up and inserts a rectangular area light source
   with size (sx, sy)
   orientation given by the normal vector (nx, ny, nz)
   centered at (tx, ty, tz)
   consisting of (lx x ly) point light sources (uniformly sampled)
   and with colour (r,g,b) - which also determines intensity

   Note that the light source is visible as a uniformly colored rectangle and
   casts no shadow. If you require a lightsource to shade another, you must
   make it into a proper solid box with backing and sides of non-light-emitting
   material
 */

  /////////////////////////////////////////////////////
  // TO DO: (Assignment 4!)
  // Implement this function to enable area light sources
  /////////////////////////////////////////////////////
  struct object3D *o;
  struct point3D *up, *forward, *right, p;
  struct pointLS *l;
  double rotation[4][4];  //rotation matrix
  double spacex, spacey;

  up = newPoint(nx, ny, nz);
  up->pw = 0;

  if (up->px == 0){
    forward = newPoint(1, 0, 0);
  }
  else{
    forward = newPoint(0, 1, 0);
  }
  forward->pw=0;

  normalize(up);
  // right = newPoint(forward->px, forward->py, forward->pz);
  // right->pw = 0;

  right = cross(up, forward);   
  free(forward);
  forward = cross(right, up);    

  rotation[0][0]=right->px;
  rotation[1][0]=right->py;
  rotation[2][0]=right->pz;
  rotation[3][0]=0;

  rotation[0][1]=forward->px;
  rotation[1][1]=forward->py;
  rotation[2][1]=forward->pz;
  rotation[3][1]=0;

  rotation[0][2]=up->px;
  rotation[1][2]=up->py;
  rotation[2][2]=up->pz;
  rotation[3][2]=0;

  rotation[0][3]=0;
  rotation[1][3]=0;
  rotation[2][3]=0;
  rotation[3][3]=1;

  free(forward);
  free(right);
  free(up);

  o=newPlane(.05,0,0,.05,1,1,1,1,1,0); // Note the plane is highly-reflective (rs=rg=.75) so we
            // ra, rd, rs, rg - Albedos for the components of the Phong model
 // r, g, b, - Colour for this plane
 // alpha - Transparency, must be set to 1 unless you are doing refraction
 // r_index - Refraction index if you are doing refraction.
 // shiny - Exponent for the specular component of the Phong model
  o->isLightSource = 1;
  Scale(o,sx,sy,1);        // Do a few transforms...
  matMult(rotation, o->T);
  Translate(o,tx+sx,ty+sy,tz);
  invert(&o->T[0][0],&o->Tinv[0][0]);    // Very important! compute
            // and store the inverse
            // transform for this object!

  insertObject(o,o_list);      // Insert into object list

  spacex = 1.0/(lx-1);
  spacey = 1.0/(ly-1);

  int j, i;

  for (j=0; j<ly; j++) {
      for (i=0; i<lx; i++) {
        p.px = i*spacex;
        p.py = j*spacey;
        p.pz = 0;
        p.pw = 1;
        matVecMult(o->T, &p);
        l = newPLS(&p, r/(lx*ly), g/(lx*ly), b/(lx*ly));
        insertPLS(l, l_list);
      }
    }
}

///////////////////////////////////
// Geometric transformation section
///////////////////////////////////

void invert(double *T, double *Tinv)
{
 // Computes the inverse of transformation matrix T.
 // the result is returned in Tinv.

 float *U, *s, *V, *rv1;
 int singFlag, i;
 float T3x3[3][3],Tinv3x3[3][3];
 double tx,ty,tz;

 // Because of the fact we're using homogeneous coordinates, we must be careful how
 // we invert the transformation matrix. What we need is the inverse of the
 // 3x3 Affine transform, and -1 * the translation component. If we just invert
 // the entire matrix, junk happens.
 // So, we need a 3x3 matrix for inversion:
 T3x3[0][0]=(float)*(T+(0*4)+0);
 T3x3[0][1]=(float)*(T+(0*4)+1);
 T3x3[0][2]=(float)*(T+(0*4)+2);
 T3x3[1][0]=(float)*(T+(1*4)+0);
 T3x3[1][1]=(float)*(T+(1*4)+1);
 T3x3[1][2]=(float)*(T+(1*4)+2);
 T3x3[2][0]=(float)*(T+(2*4)+0);
 T3x3[2][1]=(float)*(T+(2*4)+1);
 T3x3[2][2]=(float)*(T+(2*4)+2);
 // Happily, we don't need to do this often.
 // Now for the translation component:
 tx=-(*(T+(0*4)+3));
 ty=-(*(T+(1*4)+3));
 tz=-(*(T+(2*4)+3));

 // Invert the affine transform
 U=NULL;
 s=NULL;
 V=NULL;
 rv1=NULL;
 singFlag=0;

 SVD(&T3x3[0][0],3,3,&U,&s,&V,&rv1);
 if (U==NULL||s==NULL||V==NULL)
 {
  fprintf(stderr,"Error: Matrix not invertible for this object, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 // Check for singular matrices...
 for (i=0;i<3;i++) if (*(s+i)<1e-9) singFlag=1;
 if (singFlag)
 {
  fprintf(stderr,"Error: Transformation matrix is singular, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 // Compute and store inverse matrix
 InvertMatrix(U,s,V,3,&Tinv3x3[0][0]);

 // Now stuff the transform into Tinv
 *(Tinv+(0*4)+0)=(double)Tinv3x3[0][0];
 *(Tinv+(0*4)+1)=(double)Tinv3x3[0][1];
 *(Tinv+(0*4)+2)=(double)Tinv3x3[0][2];
 *(Tinv+(1*4)+0)=(double)Tinv3x3[1][0];
 *(Tinv+(1*4)+1)=(double)Tinv3x3[1][1];
 *(Tinv+(1*4)+2)=(double)Tinv3x3[1][2];
 *(Tinv+(2*4)+0)=(double)Tinv3x3[2][0];
 *(Tinv+(2*4)+1)=(double)Tinv3x3[2][1];
 *(Tinv+(2*4)+2)=(double)Tinv3x3[2][2];
 *(Tinv+(0*4)+3)=Tinv3x3[0][0]*tx + Tinv3x3[0][1]*ty + Tinv3x3[0][2]*tz;
 *(Tinv+(1*4)+3)=Tinv3x3[1][0]*tx + Tinv3x3[1][1]*ty + Tinv3x3[1][2]*tz;
 *(Tinv+(2*4)+3)=Tinv3x3[2][0]*tx + Tinv3x3[2][1]*ty + Tinv3x3[2][2]*tz;
 *(Tinv+(3*4)+3)=1;

 free(U);
 free(s);
 free(V);
}

void RotateX(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // X axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=1.0;
 R[1][1]=cos(theta);
 R[1][2]=-sin(theta);
 R[2][1]=sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateY(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Y axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][2]=sin(theta);
 R[1][1]=1.0;
 R[2][0]=-sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateZ(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Z axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][1]=-sin(theta);
 R[1][0]=sin(theta);
 R[1][1]=cos(theta);
 R[2][2]=1.0;
 R[3][3]=1.0;

 matMult(R,o->T);
}

void Translate(struct object3D *o, double tx, double ty, double tz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that translates the object by the specified amounts.

 double tr[4][4];
 memset(&tr[0][0],0,16*sizeof(double));

 tr[0][0]=1.0;
 tr[1][1]=1.0;
 tr[2][2]=1.0;
 tr[0][3]=tx;
 tr[1][3]=ty;
 tr[2][3]=tz;
 tr[3][3]=1.0;

 matMult(tr,o->T);
}

void Scale(struct object3D *o, double sx, double sy, double sz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that scales the object as indicated.

 double S[4][4];
 memset(&S[0][0],0,16*sizeof(double));

 S[0][0]=sx;
 S[1][1]=sy;
 S[2][2]=sz;
 S[3][3]=1.0;

 matMult(S,o->T);
}

void printmatrix(double mat[4][4])
{
 fprintf(stderr,"Matrix contains:\n");
 fprintf(stderr,"%f %f %f %f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[3][0],mat[3][1],mat[3][2],mat[3][3]);
}

/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize)
{
 /*
   This function sets up the camera axes and viewing direction as discussed in the
   lecture notes.
   e - Camera center
   g - Gaze direction
   up - Up vector
   fov - Fild of view in degrees
   f - focal length
 */
 struct view *c;
 struct point3D *u, *v;

 u=v=NULL;

 // Allocate space for the camera structure
 c=(struct view *)calloc(1,sizeof(struct view));
 if (c==NULL)
 {
  fprintf(stderr,"Out of memory setting up camera model!\n");
  return(NULL);
 }

 // Set up camera center and axes
 c->e.px=e->px;		// Copy camera center location, note we must make sure
 c->e.py=e->py;		// the camera center provided to this function has pw=1
 c->e.pz=e->pz;
 c->e.pw=1;

 // Set up w vector (camera's Z axis). w=-g/||g||
 c->w.px=-g->px;
 c->w.py=-g->py;
 c->w.pz=-g->pz;
 c->w.pw=0;
 normalize(&c->w);

 // Set up the horizontal direction, which must be perpenticular to w and up
 u=cross(&c->w, up);
 normalize(u);
 c->u.px=u->px;
 c->u.py=u->py;
 c->u.pz=u->pz;
 c->u.pw=0;

 // Set up the remaining direction, v=(u x w)  - Mind the signs
 v=cross(&c->u, &c->w);
 normalize(v);
 c->v.px=v->px;
 c->v.py=v->py;
 c->v.pz=v->pz;
 c->v.pw=0;

 // Copy focal length and window size parameters
 c->f=f;
 c->wl=wl;
 c->wt=wt;
 c->wsize=wsize;

 // Set up coordinate conversion matrices
 // Camera2World matrix (M_cw in the notes)
 // Mind the indexing convention [row][col]
 c->C2W[0][0]=c->u.px;
 c->C2W[1][0]=c->u.py;
 c->C2W[2][0]=c->u.pz;
 c->C2W[3][0]=0;

 c->C2W[0][1]=c->v.px;
 c->C2W[1][1]=c->v.py;
 c->C2W[2][1]=c->v.pz;
 c->C2W[3][1]=0;

 c->C2W[0][2]=c->w.px;
 c->C2W[1][2]=c->w.py;
 c->C2W[2][2]=c->w.pz;
 c->C2W[3][2]=0;

 c->C2W[0][3]=c->e.px;
 c->C2W[1][3]=c->e.py;
 c->C2W[2][3]=c->e.pz;
 c->C2W[3][3]=1;

 // World2Camera matrix (M_wc in the notes)
 // Mind the indexing convention [row][col]
 c->W2C[0][0]=c->u.px;
 c->W2C[1][0]=c->v.px;
 c->W2C[2][0]=c->w.px;
 c->W2C[3][0]=0;

 c->W2C[0][1]=c->u.py;
 c->W2C[1][1]=c->v.py;
 c->W2C[2][1]=c->w.py;
 c->W2C[3][1]=0;

 c->W2C[0][2]=c->u.pz;
 c->W2C[1][2]=c->v.pz;
 c->W2C[2][2]=c->w.pz;
 c->W2C[3][2]=0;

 c->W2C[0][3]=-dot(&c->u,&c->e);
 c->W2C[1][3]=-dot(&c->v,&c->e);
 c->W2C[2][3]=-dot(&c->w,&c->e);
 c->W2C[3][3]=1;

 free(u);
 free(v);
 return(c);
}

/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename)
{
 // Reads an image from a .ppm file. A .ppm file is a very simple image representation
 // format with a text header followed by the binary RGB data at 24bits per pixel.
 // The header has the following form:
 //
 // P6
 // # One or more comment lines preceded by '#'
 // 340 200
 // 255
 //
 // The first line 'P6' is the .ppm format identifier, this is followed by one or more
 // lines with comments, typically used to inidicate which program generated the
 // .ppm file.
 // After the comments, a line with two integer values specifies the image resolution
 // as number of pixels in x and number of pixels in y.
 // The final line of the header stores the maximum value for pixels in the image,
 // usually 255.
 // After this last header line, binary data stores the RGB values for each pixel
 // in row-major order. Each pixel requires 3 bytes ordered R, G, and B.
 //
 // NOTE: Windows file handling is rather crotchetty. You may have to change the
 //       way this file is accessed if the images are being corrupted on read
 //       on Windows.
 //
 // readPPMdata converts the image colour information to floating point. This is so that
 // the texture mapping function doesn't have to do the conversion every time
 // it is asked to return the colour at a specific location.
 //

 FILE *f;
 struct image *im;
 char line[1024];
 int sizx,sizy;
 int i;
 unsigned char *tmp;
 double *fRGB;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  f=fopen(filename,"rb+");
  if (f==NULL)
  {
   fprintf(stderr,"Unable to open file %s for reading, please check name and path\n",filename);
   free(im);
   return(NULL);
  }
  fgets(&line[0],1000,f);
  if (strcmp(&line[0],"P6\n")!=0)
  {
   fprintf(stderr,"Wrong file format, not a .ppm file or header end-of-line characters missing\n");
   free(im);
   fclose(f);
   return(NULL);
  }
  fprintf(stderr,"%s\n",line);
  // Skip over comments
  fgets(&line[0],511,f);
  while (line[0]=='#')
  {
   fprintf(stderr,"%s",line);
   fgets(&line[0],511,f);
  }
  sscanf(&line[0],"%d %d\n",&sizx,&sizy);           // Read file size
  fprintf(stderr,"nx=%d, ny=%d\n\n",sizx,sizy);
  im->sx=sizx;
  im->sy=sizy;

  fgets(&line[0],9,f);  	                // Read the remaining header line
  fprintf(stderr,"%s\n",line);
  tmp=(unsigned char *)calloc(sizx*sizy*3,sizeof(unsigned char));
  fRGB=(double *)calloc(sizx*sizy*3,sizeof(double));
  if (tmp==NULL||fRGB==NULL)
  {
   fprintf(stderr,"Out of memory allocating space for teximage\n");
   free(im);
   fclose(f);
   return(NULL);
  }

  fread(tmp,sizx*sizy*3*sizeof(unsigned char),1,f);
  fclose(f);

  // Conversion to floating point
  for (i=0; i<sizx*sizy*3; i++) *(fRGB+i)=((double)*(tmp+i))/255.0;
  free(tmp);
  im->rgbdata=(void *)fRGB;

  return(im);
 }

 fprintf(stderr,"Unable to allocate memory for image structure\n");
 return(NULL);
}

struct image *newImage(int size_x, int size_y)
{
 // Allocates and returns a new image with all zeros. Assumes 24 bit per pixel,
 // unsigned char array.
 struct image *im;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  im->sx=size_x;
  im->sy=size_y;
  im->rgbdata=(void *)calloc(size_x*size_y*3,sizeof(unsigned char));
  if (im->rgbdata!=NULL) return(im);
 }
 fprintf(stderr,"Unable to allocate memory for new image\n");
 return(NULL);
}

void imageOutput(struct image *im, const char *filename)
{
 // Writes out a .ppm file from the image data contained in 'im'.
 // Note that Windows typically doesn't know how to open .ppm
 // images. Use Gimp or any other seious image processing
 // software to display .ppm images.
 // Also, note that because of Windows file format management,
 // you may have to modify this file to get image output on
 // Windows machines to work properly.
 //
 // Assumes a 24 bit per pixel image stored as unsigned chars
 //

 FILE *f;

 if (im!=NULL)
  if (im->rgbdata!=NULL)
  {
   f=fopen(filename,"wb+");
   if (f==NULL)
   {
    fprintf(stderr,"Unable to open file %s for output! No image written\n",filename);
    return;
   }
   fprintf(f,"P6\n");
   fprintf(f,"# Output from RayTracer.c\n");
   fprintf(f,"%d %d\n",im->sx,im->sy);
   fprintf(f,"255\n");
   fwrite((unsigned char *)im->rgbdata,im->sx*im->sy*3*sizeof(unsigned char),1,f);
   fclose(f);
   return;
  }
 fprintf(stderr,"imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im)
{
 // De-allocates memory reserved for the image stored in 'im'
 if (im!=NULL)
 {
  if (im->rgbdata!=NULL) free(im->rgbdata);
  free(im);
 }
}

void cleanup(struct object3D *o_list, struct pointLS *l_list)
{
 // De-allocates memory reserved for the object list and the point light source
 // list. Note that *YOU* must de-allocate any memory reserved for images
 // rendered by the raytracer.
 struct object3D *p, *q;
 struct pointLS *r, *s;

 p=o_list;		// De-allocate all memory from objects in the list
 while(p!=NULL)
 {
  q=p->next;
  if (p->texImg!=NULL)
  {
   if (p->texImg->rgbdata!=NULL) free(p->texImg->rgbdata);
   free(p->texImg);
  }
  free(p);
  p=q;
 }

 r=l_list;
 while(r!=NULL)
 {
  s=r->next;
  free(r);
  r=s;
 }
}
