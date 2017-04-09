/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

#include "utils.h"

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
int MAX_DEPTH;

//debug globals
int intersect_count = 0;

void buildScene(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

	struct object3D *o;
	struct pointLS *l;
	struct point3D p;

 ///////////////////////////////////////
 // TO DO: For Assignment 3 you have to use
 //        the simple scene provided
 //        here, but for Assignment 4 you
 //        *MUST* define your own scene.
 //        Part of your mark will depend
 //        on how nice a scene you
 //        create. Use the simple scene
 //        provided as a sample of how to
 //        define and position objects.
 ///////////////////////////////////////

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 // Let's add a plane
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);	// Note the plane is highly-reflective (rs=rg=.75) so we
						// should see some reflections if all is done properly.
						// Colour is close to cyan, and currently the plane is
						// completely opaque (alpha=1). The refraction index is
						// meaningless since alpha=1
 Scale(o,6,6,1);				// Do a few transforms...
 RotateZ(o,PI/1.20);
 RotateX(o,PI/2.25);
 Translate(o,0,-3,10);
 invert(&o->T[0][0],&o->Tinv[0][0]);		// Very important! compute
						// and store the inverse
						// transform for this object!
 insertObject(o,&object_list);			// Insert into object list

 
 // Let's add a couple spheres

 o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6);
 Scale(o,.75,.5,1.5);
 RotateY(o,PI/2);
 Translate(o,-1.45,1.1,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,6);
 Scale(o,.5,2.0,1.0);
 RotateZ(o,PI/1.5);
 Translate(o,1.75,1.25,5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // Insert a single point light source.
 // p.px=0;
 // p.py=15.5;
 // p.pz=-5.5;
 // p.pw=1;
 // l=newPLS(&p,.95,.95,.95);
 // insertPLS(l,&light_list);
addAreaLight(3, 3, 0, 0, 5.5,\
                  0, 15.5, -5.5, 3, 3,\
                  255, 255, 255, &object_list, &light_list);


 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //

 struct colourRGB tmp_col;	// Accumulator for colour components
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj->texImg==NULL)		// Not textured, use object colour
 {
 	R=obj->col.R;
 	G=obj->col.G;
 	B=obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
 	obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }
 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////
 //col->R = obj->col.R*255;
 //col->G = obj->col.G*255;
 //col->B = obj->col.B*255;


 double NL, LN, VR; // variables used for phong
 struct point3D *L, *N, *V, *r;
	 
 struct pointLS *light = light_list;

 while(light){

	 struct point3D *light_dir = newPoint(p->px - light->p0.px, p->py - light->p0.py, p->pz - light->p0.pz); // calculating and assinging light ray
	 struct point3D *shadow_dir = newPoint(light->p0.px - p->px, light->p0.py - p->py, light->p0.pz - p->pz); // calculating and assinging light ray

	 
	 light_dir->pw = 0;
	 shadow_dir->pw = 0;

	 struct ray3D *light_ray = newRay(&(light->p0), light_dir);
	 struct ray3D *shadow_ray = newRay(p, shadow_dir);
	 
	 struct object3D *light_obj = NULL;
	 double tlambda = -1; // temp variables
	 struct point3D tp;
	 struct point3D tn;
	 double ta = 0;
	 //printf("returned: %f, %f, %f\n", light->p0.px, light->p0.py, light->p0.pz);
	 //printf("returned: %f, %f, %f\n", light_dir->px, light_dir->py, light_dir->pz);

	 findFirstHit(shadow_ray, &tlambda, obj, &light_obj, &tp, &tn, &ta, &b);
	 free(shadow_dir);
	 free(shadow_ray);
	 //printf("returned: %f, %f, %f, %f\n", p->px, p->py, p->pz, tlambda);
	 if(tlambda > 0 && tlambda < 1){
	 
	 	tmp_col.R = 0;
	 	tmp_col.G = 0;
	 	tmp_col.B = 0;

	 }


	 //phong
	 L = newPoint(light_ray->d.px, light_ray->d.py, light_ray->d.pz);
	 N = newPoint(n->px, n->py, n->pz);
	 V = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);

	 free(light_dir);
	 free(light_ray);
	 L->pw = 0;
	 N->pw = 0;
	 V->pw= 0;
	 normalize(L);
	 //normalize(N);
	 normalize(V);

	 LN = dot(L, N);

	 r = newPoint(2*LN*N->px - L->px, 2*LN*N->py - L->py, 2*LN*N->pz - L->pz);
	 r->pw = 0;
	 
	 normalize(r);

	 NL = max(0, dot(N, L));
	 //printf("NL : %f\n", L->py);
	 VR = max(0, pow(dot(V, r), obj->shinyness));
	 //printf("VR %f\n", dot(V, r));

	 double amb = obj->alb.ra;
	 double diff = obj->alb.rd * NL;
	 double spec = obj->alb.rs * VR;
	 
	 //Only multiply obj color by amb and diff, raytrace tut part2 slide 5.
	 tmp_col.R += (R*(amb + diff) + spec)*light->col.R;
	 tmp_col.G += (G*(amb + diff) + spec)*light->col.G;
	 tmp_col.B += (B*(amb + diff) + spec)*light->col.B;

	 // make sure that the colors are bound by [0, 255]

	 light = light->next;
	}

	 col->R = tmp_col.R;
	 col->G = tmp_col.G;
	 col->B = tmp_col.B;

	//printf("Depth : %i\n", depth);


	 if (depth < MAX_DEPTH){
	 	struct point3D *reflect_p = p;

	 	double dn = dot(&(ray->d), n);

	 	struct point3D *reflect_d = newPoint(ray->d.px - 2*dn*n->px, ray->d.py - 2*dn*n->py, ray->d.pz - 2*dn*n->pz);
	 	reflect_d->pw = 0;
	    //printf("returned: %f, %f, %f, %f\n", reflect_d->px, reflect_d->py, reflect_d->pz, reflect_d->pw);
	    //normalize(reflect_d);
	 	struct ray3D * reflected = newRay(reflect_p, reflect_d);
	 	rayTrace(reflected, depth+1, col, obj);

	    //free(reflect_p);
	 	free(reflect_d);
	 	free(reflected);

	 }

	 col->R = max(0, -1*max(-255, -col->R));
	 col->G = max(0, -1*max(-255, -col->G));
	 col->B = max(0, -1*max(-255, -col->B));

     free(L);
	 free(N);
	 free(V);
	 free(r);
 //fprintf(stderr,"r g b:  %f  %f\n",obj->col.R,obj->col.G);

 // Be sure to update 'col' with the final colour computed here!
	return;

}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////
	struct object3D *walker;
	double tlambda = -1; // temp variables
	struct point3D tp;
	struct point3D tn;
        //n = newPoint(1,0,0);
        //p = newPoint(0,0,0);
	if(object_list == NULL){
		return;
	}
	walker = object_list;
	while(walker != NULL){
		if(walker != Os){
			walker->intersect(walker, ray, &tlambda, &tp, &tn, a, b);
 		// 	if(*a == 3 && tlambda >= 0){
			// 	printf("returned: %f\n", tlambda);
			// } 
			if(tlambda >= 0 && (*lambda == -1 || tlambda < *lambda)){
				intersect_count += 1;
				*obj = walker;
				*lambda = tlambda;
				*p = tp;
				*n = tn;
				//printf("returned: %f, %f, %f, %f\n", tn.px, tn.py, tn.pz, tlambda);
				//n = newPoint(tn.px, tn.py, tn.pz);
                       		//printf("returned: %f, %f, %f, %f\n", p->px, p->py, p->pz, tlambda);
			}
		}
		walker = walker->next;
	}
	free(walker);
	return;
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //
 double lambda = -1;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj = NULL;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {
 	col->R=-1;
 	col->G=-1;
 	col->B=-1;
 	return;
 }

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////
 findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
 if(obj != NULL){
  //intersect_count += 1;
 	rtShade(obj, &p, &n, ray, depth, a, b, col);
 }

 //free(obj);

}

int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for raytraced pixels

  struct colourRGB background;   // Background colour
 int i,j,k,l;			// Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc<5)
 {
 	fprintf(stderr,"RayTracer: Can not parse input parameters\n");
 	fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
 	fprintf(stderr,"   size = Image size (both along x and y)\n");
 	fprintf(stderr,"   rec_depth = Recursion depth\n");
 	fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
 	fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
 	exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
 	fprintf(stderr,"Unable to allocate memory for raytraced image\n");
 	exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 3, you can use
 //        the simple scene already provided. But
 //        for Assignment 4 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();		// Create a scene. This defines all the
			// objects in the world of the raytracer

 //////////////////////////////////////////
 // TO DO: For Assignment 3 you can use the setup
 //        already provided here. For Assignment 4
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-3;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
 g.px=0;
 g.py=0;
 g.pz=3;
 g.pw=0;

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=0;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
 	fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
 	cleanup(object_list,light_list);
 	deleteImage(im);
 	exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 struct colourRGB supersampledColor;
 struct point3D * imagePlane;
 struct point3D * origin;
 struct point3D * direction;
 //set to 1 if you want no anti-aliasing
 int supersamplingSize = 4;
 if(!antialiasing) supersamplingSize = 1;
 fprintf(stderr,"Rendering row: ");
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
 	fprintf(stderr,"%d/%d, ",j,sx);
 	for (i=0;i<sx;i++)
 	{
    ///////////////////////////////////////////////////////////////////
    // TO DO - complete the code that should be in this loop to do the
    //         raytracing!
    ///////////////////////////////////////////////////////////////////
		origin  = newPoint(0,0,0);//newPoint(cam->e.px,cam->e.py,cam->e.pz);
		imagePlane = newPoint(0,0,0);
 		/*pc.px = 0;
		pc.py = 0;
		pc.pz = 0;
		pc.pw = 0;
		*/
		col.R = 0; col.G = 0; col.B = 0;
		supersampledColor.R = 0; supersampledColor.G = 0; supersampledColor.B = 0;
		l=0;k=0;	
		for (l=0;l<supersamplingSize;l++)
		{
			for(k=0;k<supersamplingSize;k++)
			{
				double dx = (0.5 + (1.0*k))/supersamplingSize;
				double dy = (0.5 + (1.0*l))/supersamplingSize;
				imagePlane->px = (cam->wl) + (i + dx)*du;
				imagePlane->py = (cam->wt) + (j + dy)*dv;
				imagePlane->pz = cam->f;
 				/*d.px = (cam->wl) + (i + dx)*du;
                                d.py = (cam->wt) + (j + dy)*dv;
                                d.pz = cam->f;
			        d.pw = 0;*/ 
		                //fprintf(stderr,"(%d, %d): %f/%f, \n",l, k, d.px, d.py);
		 		//Ray Direction 
				direction = newPoint(imagePlane->px, imagePlane->py, imagePlane->pz); 	
				subVectors(origin, direction);
				//fprintf(stderr,"(%d, %d): %f/%f, \n",1, 1, d.px, d.py);
				direction->pw = 0;
			    	//Convert to world-space
				matVecMult(cam->C2W, direction);
				matVecMult(cam->C2W, origin);
				//matVecMult(cam->C2W, &d);
                                //matVecMult(cam->C2W, &pc);
				direction->pw = 0;
				d.pw = 0;
				ray = newRay(origin, direction);
				//ray = newRay(&pc, &d);
				rayTrace(ray, 0, &col, NULL);
				//weight the color using a gaussian
				//double colorWeightk = (1/(0.5*sqrt(2*PI))) * exp(-1*(pow((k-((supersamplingSize-1)/2)),2))/(2*0.25));
				//double colorWeightl = (1/(0.5*sqrt(2*PI))) * exp(-1*(pow((l-((supersamplingSize-1)/2)),2))/(2*0.25));
				double colorWeight = 1;//(-1*(max(-1*colorWeightk, -1*colorWeightl)));
				//printf("d: %f, p:%f, result: %f\n",pow(k,2), -1*(pow((k-(supersamplingSize/2)),2))/(2*0.25), colorWeightk);
				supersampledColor.R += col.R*colorWeight;
				supersampledColor.G += col.G*colorWeight;
				supersampledColor.B += col.B*colorWeight; 
				free(direction);
			}
		}
		col.R = supersampledColor.R/(supersamplingSize*supersamplingSize);
		col.G = supersampledColor.G/(supersamplingSize*supersamplingSize);
		col.B = supersampledColor.B/(supersamplingSize*supersamplingSize);
		//printf("%f, %f, %f\n", col.R, col.G, col.B);
		 //Paint the RGB
		*rgbIm  = col.R;
		rgbIm++;
		*rgbIm  = col.G;
		rgbIm++;
		*rgbIm  = col.B;
		rgbIm++;
		free(origin);
 		free(imagePlane);

  	} // end for i
 } // end for j
 
 
 

 fprintf(stderr,"%d, ", intersect_count);
 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list);		// Object and light lists
 deleteImage(im);				// Rendered image
 free(cam);					// camera view
 exit(0);
}

