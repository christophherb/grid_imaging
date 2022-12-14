/**
\mainpage
Simple Meta-Conic Neutron Raytracer is a framework for raytracing geometries of the form: @f$ r^2=k_1 + k_2 z + k_3 z^2 @f$.

<h3>General Notes</h3>
To use the software you must make a Scene element using the function makeScene(). You must then add items to this scene element using the various add function (addDisk(), addParaboloid(), etc...). Next you must call the function traceSingleNeutron() for every neutron you would like to trace through the geometry. The maximum number of each geometry you can place in a scene is defined by the MAX_CONICSURF, MAX_DISK and MAX_DETECTOR definitions in the conic.h file.

<h3>TODO</h3>

@todo
       Name variable for each component <br/>
       Normalize Detector Events by weight of neutron <br/>

<h3>Known Bugs</h3>
@bug  HPPH works incorrectly for M != 1 <br/>
      Neutrons t=1e-11 away from a surface will pass through <br/>

<h3>Note on Pointers</h3>
This framework uses pointers extensivly. It is important to be familiar with how to use pointers.
<br/>Here is a quick overview.

@code
//Making an Element
ConicSurf c = makeParaboloid(...);
int k = 10;

//Making a Pointer from an Element
ConicSurf * c_pointer = &c;
int * k_pointer = &k;

//Making an Element from a Pointer
ConicSurf c2 = *c_pointer;
int ten = *k_pointer;

//Getting Item in Element
double k1 = c.k1;

//Getting Item in Element from a Pointer
double k1 = c_pointer->k1;

//Functions that have pointers as parameters can modify the element being pointed to.
//Functions that have elements as parameters can not modify the value of the element.
@endcode

<h3>Stand Alone Example Code</h3>
This framework can be used entirely by itself as the following example demonstrates.
@code
/////////////////////////////////////////////////////////////////////////////////
// Giacomo Resta <gresta@mit.edu>
//
// Basic standalone exmaple of a raytracer. It is advised to modify
// the getRandom() function to use a better random number generator (i.e. glib).
// The systems random generator was used to preserve clarity and conciseness.
/////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>

#include "conic.h"
#include "w1_general.h"

#define NUM_NEUTRON 1000000
*/
//Function to get random number
double getRandom() {
    return rand01();
    //return (double)lrand48()/RAND_MAX;
}
/*
//Function to get new particle from source
Particle generateParticleFromSource(double radius, double div, double Lambda0,
        double num_neutrons) {

    double chi = 2*M_PI*(getRandom());
    double r = sqrt(getRandom()) * radius;

    double v = 3956.036 / Lambda0;

    double theta = sqrt(getRandom())*div;
    double phi = getRandom()*2*M_PI;
    double tan_r = tan(theta);

    double vz = v/sqrt(1+tan_r*tan_r);
    double vr = tan_r * vz;

    return makeParticle(r*cos(chi),r*sin(chi),0,cos(phi)*vr,sin(phi)*vr,
                vz,0,0.0,0.0,0.0,1.0/num_neutrons);
}

//Function to add items to scene
void addItems(Scene* s, double instr_len,double r,double f,double M,
        double max_mir_len,double m, double mirr_thick) {

    //Change code here to make different Scenes
    PP p = addPPShell(0.0, instr_len, r, f, M, max_mir_len, m, m, mirr_thick, s);
    addDisk(p.p0->zs, 0.0, rConic(p.p0->ze, *p.p0)*p.p0->zs/p.p0->ze, s);
    addDisk(p.p1->zs, rConic(p.p1->zs, *p.p1), 10000,s);
    addDetector(10.0, -0.01, 0.01, -0.01, 0.01, 600, 600, NUM_NEUTRON, "test.txt", s);
}

//Main Function
int main() {
    //seed random generator (starts the generator)
    srand48((long)time(0));

    //Make a new Scene
    Scene s = makeScene();

    //Add Items and initialize Scene
    addItems(&s,10.0,0.068,4.2,1,0.7,3,0.001);
    initSimulation(&s);

    //Raytrace all particles through the Scene
    double i;
    for (i = 0; i < NUM_NEUTRON; i++) {
        Particle p = generateParticleFromSource(0.005, 0.02422, 4, NUM_NEUTRON);
        traceSingleNeutron(&p,s);
    }

    //Finish Simulation of the Scene
    finishSimulation(&s);

    return 0;
}
@endcode

*/

/**
    @file conic.h
    \brief General Framework for generating and raytracing geometries of the
    form @f$ r = k_1+k_2 z+k_3 z^2 @f$

    @author Giacomo Resta <gresta@mit.edu>
    @version 0.2

    @section LICENSE
    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the Software
    is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
    PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
    FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

    @section DESCRIPTION
     General Framework for generating and raytracing geometries of the form
     @f$ r = k_1+k_2 z+k_3 z^2 @f$
*/

/////////////////////////////////////
// Simulation
/////////////////////////////////////

#ifndef MIT_CONICS
#define MIT_CONICS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** @defgroup simgroup Simulator Internals
    Contains items general to the simulation
    @{
*/
//! Max number of ConicSurf allowed in a Scene
#define MAX_FLATSURF 200

//! Max number of ConicSurf allowed in a Scene
#define MAX_CONICSURF 100

//! Max number of Disks allowed in a Scene
#define MAX_DISK 100

//! Max number of Detectors allowed in a Scene
#define MAX_DETECTOR 10

//! If "1" simulator will record z location where neutron with greatest grazing angle reflected for each ConicSurf
/*! The information is stored in the max_ga and max_ga_z0 members of each ConicSurf, which are only present if
    this flag is 1. See source code for clarification.

@note You must use default traceNeutronConic function or implement saving routine yourself */
#define REC_MAX_GA 0

#define V2Q_conic 1.58825361e-3
#define Q2V_conic 629.622368
#define POT_V 46.839498800356665//theta_crit = 4 pi^2 kappa^2/m_n^2
#define m_Si 0.478 // m-value of pure silicon
//! Stucture to represent a point
typedef struct {
    double x; //!< x-axis position of point
    double y; //!< y-axis position of point
    double z; //!< z-axis position of point
} Point;

//! Structure to represent a vector
typedef struct {
    double x; //!< x-axis length of vector
    double y; //!< y-axis length of vector
    double z; //!< z-axis length of vector
} Vec;

//! Structure to represent a particle
typedef struct {
    double _x; //!< x axis position of particle
    double _y; //!< y axis position of particle
    double _z; //!< z axis position of particle
    double _vx; //!< x axis components of velocity
    double _vy; //!< y axis components of velocity
    double _vz; //!< z axis components of velocity
    double _sx; //!< x spin vector components
    double _sy; //!< y spin vector components
    double _sz; //!< z spin vector components
    double w; //!< Weight of particle
    int silicon; //!< +1 if Particle is in silicon -1 if Particle is in air
    int absorb; //!< Absorb flag (0 is not absorbed)
    double _t; //!< Time of flight of particle
} Particle;

/*! \brief Function to make a point

@param x x-axis position
@param y y-axis position
@param z z-axis position
*/
Point makePoint(double x, double y, double z) {
    Point p;
    p.x = x;
    p.y = y;
    p.z = z;

    return p;
}

/*! \brief Function to make a vector

@param x x-axis length
@param y y-axis length
@param z z-axis length
*/
Vec makeVec(double x, double y, double z) {
    Vec p;
    p.x = x;
    p.y = y;
    p.z = z;

    return p;
}

//! Function to compute length of a vector
double getMagVec(Vec v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

//! Function to compute dot product of two vectors
double dotVec(Vec v1, Vec v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//! Function to compute the sum of two vectors
Vec difVec(Vec v1, Vec v2){
    return makeVec(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}

//! Function to compute the sum of two vectors
Vec sumVec(Vec v1, Vec v2){
    return makeVec(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}

//! Function to compute the sum of two vectors
Vec skalarVec(Vec v1, double k){
    return makeVec(v1.x*k, v1.y*k, v1.z*k);
}
/*! \brief Function to make a particle

@param x x-axis position
@param y y-axis position
@param z z-axis position
@param vx x-axis velocity
@param vy y-axis velocity
@param vz z-axis velocity
@param t time of flight of neutron
@param sx x-axis component of spin vector
@param sy y-axis component of spin vector
@param sz z-axis component of spin vector
@param w weight of particle
*/
Particle makeParticle(double x, double y, double z,
    double vx, double vy, double vz, double t,
    double sx, double sy, double sz, int silicon, double w) {

    Particle pa;

    pa._x = x;
    pa._y = y;
    pa._z = z;

    pa._vx = vx;
    pa._vy = vy;
    pa._vz = vz;

    pa._sx = sx;
    pa._sy = sy;
    pa._sz = sz;
    pa.silicon = silicon;
    pa.w = w;

    pa.absorb = 0;
    pa._t = t;

    return pa;
}

//! Function to get position of particle
Point getParticlePos(Particle p) {
    return makePoint(p._x,p._y,p._z);
}

//! Function to get velocity vector of particle
Vec getParticleVel(Particle p) {
    return makeVec(p._vx, p._vy, p._vz);
}

/*! \brief Function to move particle a specific time step.

Will not move particle if t < 0. Does not simulate
gravity.

@param t time step to move particle
@param p pointer to particle
*/
void moveParticleT(double t, Particle* p) {
    if (t < 0)
        return;
    if (p->silicon>0){
        p->_x = p->_x+p->_vx*t;
        p->_y = p->_y+p->_vy*t;
        p->_z = p->_z+p->_vz*t;
        p->_t = p->_t+t;
        p->w *= exp(-t*98900/52.338);//ugly hard coding in the penetration depth of silicon at 989 m/s
        //absorbParticle(p);
    }
    else{
        p->_x = p->_x+p->_vx*t;
        p->_y = p->_y+p->_vy*t;
        p->_z = p->_z+p->_vz*t;
        p->_t = p->_t+t;
    }
}

/*! \brief Function to move particle to position z.

Will not move particle if moving particle to z
position would require negative time.
Does not simulate gravity.

@param z z-position to move particle
@param p pointer to particle
*/
void moveParticleZ(double z, Particle* p) {
    double t = (z-p->_z)/p->_vz;
    moveParticleT(t, p);
}

/*! \brief Function to compute new position of particle
without modifying the position of the actual particle.

Will not move particle if t < 0. Does not simulate gravity.

@param t timestep to move particle
@param p particle
*/
Particle copyMoveParticleT(double t, Particle p) {
    Particle p2 = p;
    moveParticleT(t,&p2);
    return p2;
}

/*! \brief Function to move particle to position z
without modifying the position of the actual particle.

Will not move particle if moving particle to z
position would require negative time.
Does not simulate gravity.

@param z z-position to move particle
@param p pointer to particle
*/
Particle copyMoveParticleZ(double z, Particle p) {
    Particle p2 = p;
    moveParticleZ(z,&p2);
    return p2;
}

/*! \brief Mathematical Aid for Snell's Law for reflection.

Does not take into account grazing angle constrictions.
Only computes mathematical reflection.

@param n Normal vector
@param p Pointer to particle
*/
void reflectParticle(Vec n, Particle* p) {
    double vn = dotVec(getParticleVel(*p),n);

    p->_vx = p->_vx-2*vn*n.x;
    p->_vy = p->_vy-2*vn*n.y;
    p->_vz = p->_vz-2*vn*n.z;
}

/*! \brief Function to mark particle as absorbed

@param p Pointer to particle to be absorbed */
void absorbParticle(Particle* p)  {
    p->_vx = 0;
    p->_vy = 0;
    p->_vz = 0;
    p->w = 0;
    p->absorb = 1;
}

/*! \brief Function to set weight of particle.

Will set the weight of the particle to w.  */
void setWeightParticle(double w, Particle* pa) {
    pa->w = w;
}

/*! \brief Function to solve quadratic equations for smallest positive result.

If no positive result returns -1. Parameters are coefficents such that
@f$ 0=A z^2 + B z + C @f$

@return Will return smallest positive value or -1 if no smallest positive value
*/
double solveQuad(double A, double B, double C) {
    if (fabs(A) < 1e-11 && B != 0)           //FIXME: 1e-11 cutoff may cause problems
        return -C/B;
    else {
        double det = B*B - 4*A*C;
        if (det < 0)
            return -1;
        else {
            double sdet = sqrt(det);
            double s1 = (-B+sdet)/(2*A);
            double s2 = (-B-sdet)/(2*A);

            if (fabs(s1) < 1e-11) s1=0.0;     //FIXME: 1e-11 cutoff may cause problems
            if (fabs(s2) < 1e-11) s2=0.0;     //FIXME: 1e-11 cutoff may cause problems

            if (s1 > 0.0) {
                if (s2 > 0.0) {
                    if (s1 > s2)
                        return s2;
                    return s1;
                }
                else
                    return s1;
           }
           if (s2 > 0.0)
               return s2;
        }
    }
    return -1;
}

//! Returns sign of x, either -1 or 1
int sign(double x) {
    if (x < 0)
        return -1;
    return 1;
}

/** @} */ //end of simgroup

/*! @ingroup detectorgroup
\brief Structure for representing inline detectors

@warning Do not directly modify this structure*/
typedef struct {
    double z0; //!< z-axis position of detector
    double xmin; //!< Smallest x value to detect
    double xmax; //!< Largest x value to detect
    double xstep; //!< x size of subsampling pixel
    double ymin; //!< Smallest y value to detect
    double ymax; //!< Largest y value to detect
    double ystep; //!< y size of subsampling pixel
    int nrows;    //!< Number of pixels along y axis
    int ncols;    //!< Number of pixels along x axis
    double num_particles; //!< Number of particles being emitted from source
    double *num_count; //!< Pointer to the number of particles that hit detector
    double **data; //!< Pointer to 2d data, pixel_x_y = data[x][y]
    char* filename; //!< Name of output file of detector (should end in .txt)
} Detector;

/*! @ingroup diskgroup
\brief Structure for representing Disk geometry

 Creates a doughnut with inner radius r0 and outer radus r1 at position z0.
Neutrons between r0 and r1 are absorbed. */
typedef struct {
    double r0; //!< Inner radius of doughnut
    double r1; //!< Outer radius of doughnut
    double z0; //!< z-axis position of Disk
} Disk;

/*! @ingroup conicgroup */
enum ConicType {
    PARA,
    HYPER,
    ELLIP
};


/*! @ingroup conicgroup
\brief Structure to contain z-axis symetric conic sections

Contains any geometry that can be expressed as
@f$ r^2=k_1 + k_2 z + k_3 z^2 @f$

@warning Do not directly modify values in this structure directly */
typedef struct {
    double k1; //!< @f$ k_1 @f$ in equation below
    double k2; //!< @f$ k_2 @f$ in equation below
    double k3; //!< @f$ k_3 @f$ in equation below
    double zs; //!< z-axis position of start of mirror
    double ze; //!< z-axis position of end of mirror
    double m;  //!< m value for mirror (1.0 for Nickel)
    int doubleReflections; //!< 0 if reflections from the back of the surface cannot happen 1 otherwise
    //Only for reference
    double f1; //!< z-axis position of first focus
    double f2; //!< z-axis position of second focus, for paraboloid this is unassigned
    double a;  //!< Value of a, specific to geometry type
    double c;  //!< Value of c, for paraboloid this is unassigned
    enum ConicType type; //!< Type of mirror geometry

    #if REC_MAX_GA
    double max_ga; //!< Max Grazing Angle of Reflected Neutron (Exists only if REC_MAX_GA)
    double max_ga_z0; //!< Collision point of Max Grazing Neutron (Exists only if REC_MAX_GA)
    #endif

} ConicSurf;


/*! @ingroup flatgroup
\brief Structure to contain z-axis symetric flat sections

Contains flat geometries which can be expressed as
@f$ x^2 = k_1 + k_2 z + k_3 z^2 @f$ or
@f$ y^2 = k_1 + k_2 z + k_3 z^2 @f$

@warning Do not directly modify values in this structure directly */
typedef struct {
    double k1; //!< @f$ k_1 @f$ in equation below
    double k2; //!< @f$ k_2 @f$ in equation below
    double k3; //!< @f$ k_3 @f$ in equation below
    double zs; //!< z-axis position of start of mirror
    double ze; //!< z-axis position of end of mirror
    double ll; //!< left/lower limit of mirror along translational symmetry
    double rl; //!< right/upper limit of mirror along translational symmetry
    double m;  //!< m value for mirror (1.0 for Nickel)

    //Only for reference
    double f1; //!< z-axis position of first focus
    double f2; //!< z-axis position of second focus, for paraboloid this is unassigned
    double a;  //!< Value of a, specific to geometry type
    double c;  //!< Value of c, for paraboloid this is unassigned
    //enum FlatType type; //!< Type of mirror geometry
    int doubleReflections; // will determine whether the geometry allows double reflections
    #if REC_MAX_GA
    double max_ga; //!< Max Grazing Angle of Reflected Neutron (Exists only if REC_MAX_GA)
    double max_ga_z0; //!< Collision point of Max Grazing Neutron (Exists only if REC_MAX_GA)
    #endif

} FlatSurf;
/*! @ingroup simgroup
\brief Structure to hold all scene geometry

The number of possible ConicSurf, Disk and Detector in the Scene are
determined by MAX_CONICSURF, MAX_DISK and MAX_DETECTOR.
*/
typedef struct {
    FlatSurf f[MAX_FLATSURF]; //!< Array of all ConicSurf in Scene
    int num_f;                  //!< Number of ConicSurf in Scene

    ConicSurf c[MAX_CONICSURF]; //!< Array of all ConicSurf in Scene
    int num_c;                  //!< Number of ConicSurf in Scene

    Disk di[MAX_DISK];          //!< Array of all Disk in Scene
    int num_di;                 //!< Number of Disk in Scene

    Detector d[MAX_DETECTOR];  //!< Array of all Detector in Scene
    int num_d;                 //!< Number of Detector in Scene

    //! Function called to handle Neutron-Flat Interaction
    void (*traceNeutronFlat)(Particle*,FlatSurf);
    //! Function called to handle Neutron-Conic Interaction
    void (*traceNeutronConic)(Particle*,ConicSurf);
    //! Function called to handle Neutron-Disk Interaction
    void (*traceNeutronDisk)(Particle*,Disk);
    //! Function called to handle Neutron-Detector Interaction
    void (*traceNeutronDetector)(Particle*,Detector);

} Scene;

/////////////////////////////////////
// Inline Detector
/////////////////////////////////////

/** @defgroup detectorgroup Detector
    Contains code related to the inline detectors
    @{
*/

/*! \brief Function to make Detector

@param z0 z-axis position of detector
@param xmin Smallest x value to detect
@param xmax Largest x value to detect
@param ymin Smallest y value to detect
@param ymax Largest y value to detect
@param xres Number of pixels along x axis
@param yres Number of pixels along y axis
@param num_particles Total number of particles being emitted
@param filename Name of output file of detector (should end in .txt)
*/
Detector makeDetector(double z0,double xmin, double xmax, double ymin, double ymax, int xres,
    int yres, double num_particles, char* filename) {

    Detector d;
    d.z0 = z0;
    d.xmin = xmin;
    d.xmax = xmax;
    d.xstep = (xmax-xmin)/xres;

    d.ymin = ymin;
    d.ymax = ymax;
    d.ystep = (ymax-ymin)/yres;

    d.ncols = xres;
    d.nrows = yres;
    d.filename = filename;
    d.num_particles = num_particles;

    d.num_count = (double*)malloc(sizeof(double));
    if (d.num_count == NULL) {
        fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
        exit(-1);
    }
    (*d.num_count) = 0;

    d.data = (double**)malloc(d.ncols*sizeof(double *));
    if (d.data == NULL) {
        fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
        exit(-1);
    }
    int x;
    for(x = 0; x  < d.ncols; x++) {
        d.data[x] = (double*)malloc(d.ncols*sizeof(double));
        if (d.data[x] == NULL) {
            fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
            exit(-1);
        }
        (*d.data[x]) = 0;
    }

    return d;
}

/*! \brief Function to make and add Detector

@param z0 z-axis position of detector
@param xmin Smallest x value to detect
@param xmax Largest x value to detect
@param ymin Smallest y value to detect
@param ymax Largest y value to detect
@param xres Number of pixels along x axis
@param yres Number of pixels along y axis
@param num_particles Total number of particles being emitted
@param filename Name of output file of detector (should end in .txt)
@param s Scene to add Detector to
*/
Detector* addDetector(double z0, double xmin, double xmax, double ymin, double ymax, double xres,
    double yres, double num_particles, char* filename, Scene* s) {
    if (s->num_d >= MAX_DETECTOR-1) {
        fprintf(stderr,"TOO MANY DETECTORS IN SCENE");
        exit(-1);
    }
    s->d[s->num_d] = makeDetector(z0,xmin,xmax,ymin,ymax,xres,yres,num_particles,filename);
    s->num_d++;
    return &s->d[s->num_d-1];
}

/*! \brief Function to compute time of first collision for a Detector.

@param p Particle to consider
@param d Detector to consider

@return Time until the propogation or -1 if particle will not hit detector
*/
double getTimeOfFirstCollisionDetector(Particle p, Detector d) {
    double t = (d.z0-p._z)/p._vz;
    if (t <= 0)
        return -1;
    Particle p2 = copyMoveParticleT(t,p);
    if (p2._x > d.xmax || p2._x < d.xmin || p2._y > d.ymax || p2._y < d.ymin)
        return -1;
    return t;
}

/*! \brief Function to raytrace Detector

@param p Pointer to particle to be traced
@param d Detector to be traced
*/
void traceNeutronDetector(Particle* p, Detector d) {
    double t = getTimeOfFirstCollisionDetector(*p, d);
    if (t < 0)
        return;
    moveParticleT(t,p);
    d.data[(int)floor((p->_x-d.xmin)/d.xstep)][(int)floor((p->_y-d.ymin)/d.ystep)] += p->w;
    (*d.num_count) += p->w;
}

/*! \brief Function to finalize detector

Will write data and free data array.

@param d Detector to finalize
*/
void finishDetector(Detector d) {
    int x,y;
    if (d.filename != "") {
        FILE *file;
        file = fopen(d.filename,"w");

        double intensity = (*d.num_count);
        fprintf(file, "#I=%e I_ERR=%e xmin=%f xmax=%f ymin=%f ymax=%f ncols=%i nrows=%i\n",
            intensity, sqrt(intensity/d.num_particles), d.xmin, d.xmax, d.ymin, d.ymax, d.ncols, d.nrows); //FIXME: check I_ERR sqrt(I/num_particles)

        //Write data
        for (x=0; x < d.ncols; x++) {
            for (y=0; y < d.nrows; y++)
                fprintf(file, "%e ", d.data[x][y]);
            fprintf(file, "\n");
        }
        fclose(file);
    }
    for (x=0; x < d.ncols; x++)
        free(d.data[x]);
    free(d.data);
    free(d.num_count);
}

/** @} */ //end of detectorgroup

/////////////////////////////////////
// Geometry Types
/////////////////////////////////////

/////////////////////////////////////
// Disks
/////////////////////////////////////

/** @defgroup diskgroup Disk
    Contains code related to Disks
    @{
*/

/*! \brief Function for creating a Disk structure

@param z0 z-axis position of Disk
@param r0 Inner radius of doughnut
@param r1 Outer radius of doughnut

@see Disk
*/
Disk makeDisk(double z0, double r0, double r1) {
    Disk d;

    d.r0 = r0;
    d.z0 = z0;
    d.r1 = r1;

    return d;
}

/*! \brief Function for making and adding Disk to Scene

@param z0 z-axis position of Disk
@param r0 Inner radius of doughnut
@param r1 Outer radius of doughnut
@param s Scene to add Disk to

@see Disk
*/
Disk* addDisk(double z0, double r0, double r1, Scene* s) {
    if (s->num_di >= MAX_DISK-1) {
        fprintf(stderr,"TOO MANY DISKS IN SCENE");
        exit(-1);
    }
    s->di[s->num_di] = makeDisk(z0, r0, r1);
    s->num_di++;
    return &s->di[s->num_di -1];
}

/*! \brief Function to compute time of first collision for a disk

@param p Particle to consider
@param d Disk to consider
@return Time until the propogation or -1 if particle will not hit disk
*/
double getTimeOfFirstCollisionDisk(Particle p, Disk d) {
    double tz = (d.z0-p._z)/p._vz;
    if (tz <= 0)
        return -1;
    Particle p2 = copyMoveParticleT(tz, p);
    double rp = sqrt(p2._x*p2._x+p2._y*p2._y);
    if (rp > d.r0 && rp < d.r1 && fabs(p2._z-d.z0) < 1e-11)
        return (d.z0-p._z)/p._vz;
    return -1;
}

/*! \brief Function to raytrace Disks

@param p Pointer to particle to be traced
@param d Disk to be traced
*/
void traceNeutronDisk(Particle* p, Disk d) {
    double t = getTimeOfFirstCollisionDisk(*p, d);

    if (t <= 0)
        return;

    moveParticleT(t, p);
    //absorbParticle(p); //Disk will only be used to propagate neutrons somewhere
}

/** @} */ //end of diskgroup

/////////////////////////////////////
// Z-Axis Symetric Conic Sections
/////////////////////////////////////

/** @defgroup conicgroup ConicSurf
    Contains code related to ConicSurfs
    @{
*/

/*! \brief Function to return radius of ConicSurf at a z-axis position.

Will return radius even if z is outside the bounds of zs and ze
for the particular ConicSurf.

@param z z-axis position to compute radius
@param s ConicSurf to compute radius of
*/
double rConic(double z, ConicSurf s) {
    return sqrt(s.k1+s.k2*z+s.k3*z*z);
}

/*! \brief Function for generating Hyperboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeHyperboloid(double f1, double f2, Point p,
   double zstart, double zend, double m) {
    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;

    double r2 = p.x*p.x+p.y*p.y;
    double c = (f1-f2)/2;

    double u = p.z+c-f1;
    double a = sqrt(((u*u+c*c+r2)-sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

    s.k3 = c*c/(a*a)-1;
    s.k2 = 2*s.k3*(c-f1);
    s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

    s.m = m;
    s.f1 = f1;
    s.f2 = f2;
    s.a = a;
    s.c = c;

    s.type = HYPER;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating Ellipsoid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeEllipsoid(double f1, double f2, Point p,
    double zstart, double zend, double m, int doubleReflections) {
    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double c = (f1-f2)/2;

    double u = p.z+c-f1;
    double a = sqrt(((u*u+c*c+r2)+sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

    s.k3 = c*c/(a*a)-1;
    s.k2 = 2*s.k3*(c-f1);
    s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

    s.m = m;
    s.f1 = f1;
    s.f2 = f2;
    s.a = a;
    s.c = c;

    s.type = ELLIP;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}


/*! \brief Function for generating Flat Ellipse with symmetry along the vertical y.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param b the short half axis of the ellipse, positive for translational symmetry along y, negative for translational symmetry along x
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the left/lower limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param rl the right/upper limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param m m value for reflectivity of the surface


@see ConicSurf
*/
FlatSurf makeFlatEllipse(double f1, double f2, Point p, double zstart, double zend, double ll, double rl, double m, int doubleReflections) {
    FlatSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    s.ll = ll;
    s.rl = rl;
    double r2 = p.x*p.x + p.y*p.y;
    double c = (f1-f2)/2;
    double u = p.z+c-f1;
    double a = sqrt(((u*u+c*c+r2)+sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

    s.k3 = c*c/(a*a)-1;
    s.k2 = 2*s.k3*(c-f1);
    s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

    s.m = m;
    s.f1 = f1;
    s.f2 = f2;
    s.a = a;
    s.c = c;

    //s.type = ELLIP;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating Paraboloid ConicSurf.

@param f z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeParaboloid(double f, Point p, double zstart,
    double zend, double m, int doubleReflections) {

    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double a = (-(p.z-f)+sign(p.z-f)*sqrt((p.z-f)*(p.z-f)+r2))/2;

    s.k3 = 0.0;
    s.k2 = 4*a;
    s.k1 = s.k2*(a-f);

    s.m = m;
    s.f1 = f;
    s.a = a;

    s.type = PARA;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating Flat Parabola for FlatSurf.

@param f z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror, putting one of x or y to 0 results in the surface being parallel to said coordinate
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the left/lower limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param rl the right/upper limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param m m value for reflectivity of the surface
@param doubleReflections wether double reflections are allowed

@see FlatSurf
*/
FlatSurf makeFlatparbola(
    double f,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections) {

    FlatSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double a = (-(p.z-f)+sign(p.z-f)*sqrt((p.z-f)*(p.z-f)+r2))/2;

    s.k3 = 0.0;
    s.k2 = 4*a;
    s.k1 = s.k2*(a-f);

    s.m = m;
    s.f1 = f;
    s.a = a;
    s.ll = ll;
    s.rl = rl;
    //s.type = PARA;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating and adding Paraboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Paraboloid to

@see ConicSurf
*/
ConicSurf* addParaboloid(double f1, Point p, double zstart, double zend,
    double m, int doubleReflections, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeParaboloid(f1,p,zstart,zend,m, doubleReflections);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding a flat Parabolic FlatSurf.

@param f1 z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the lower bound of the mirror
@param rl the upper bound of the mirror
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to
@param doubleReflections wether double reflections can occur
@see FlatSurf
*/
FlatSurf* addFlatParabola(
    double f,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    Scene* s,
    int doubleReflections) {
    if (s->num_f >= MAX_FLATSURF-1) {
        fprintf(stderr,"TOO MANY FLATSURF IN SCENE");
        exit(-1);
    }
    s->f[s->num_f] = makeFlatparbola(f,p,zstart,zend,ll,rl,m,doubleReflections);
    s->num_f++;
    return &s->f[s->num_f-1];
}

/*! \brief Function for generating and adding Hyperboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Hyperboloid to

@see ConicSurf
*/
ConicSurf* addHyperboloid(double f1, double f2, Point p, double zstart,
    double zend, double m, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeHyperboloid(f1,f2,p,zstart,zend,m);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding Ellipsoid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to

@see ConicSurf
*/
ConicSurf* addEllipsoid(double f1, double f2, Point p, double zstart,
    double zend, double m, int doubleReflections, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeEllipsoid(f1,f2,p,zstart,zend,m,doubleReflections);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding a flat Ellipse FlatSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to

@see ConicSurf
*/
FlatSurf* addFlatEllipse(
    double f1,
    double f2,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections,
    Scene* s) {
    if (s->num_f >= MAX_FLATSURF-1) {
        fprintf(stderr,"TOO MANY FLATSURF IN SCENE");
        exit(-1);
    }
    s->f[s->num_f] = makeFlatEllipse(f1,f2,p,zstart,zend,ll,rl,m,doubleReflections);
    s->num_f++;
    return &s->f[s->num_f-1];
}
//!TODO
double getGrazeAngleConic(Particle p, ConicSurf s) {
    /*
    double v = sqrt(dotVec(getParticleVel(p),getParticleVel(p)));
    double vn = dotVec(getParticleVel(p),n);
    return fabs(acos(vn/v)) - M_PI/2;
    */
}

/*! \brief Function for returning normal vector of ConicSurf at Point p

Will compute vector even if p is not on surface.
MAKE SURE p IS ON SURFACE

@param p Point to compute normal vector
@param s ConicSurf to compute normal vector of
*/
Vec getNormConic(Point p, ConicSurf s) {
    double det = s.k2*s.k2+4*s.k3*(p.x*p.x+p.y*p.y-s.k1);
    if (det <= 0.){

        return makeVec(-p.x/sqrt(p.x*p.x + p.y*p.y),-p.y/(p.x*p.x + p.y*p.y),0);
    }
    double den = sqrt(det);
    double nx = -2*p.x/den;
    double ny = -2*p.y/den;
    double nz = sign(2*s.k3*p.z+s.k2);
    double n = sqrt(nx*nx+ny*ny+nz*nz);
    //printf("%f,%f,%f \n", nx/n,ny/n,nz/n);
    return makeVec(nx/n,ny/n,nz/n);
}

/*! \brief Function for returning normal vector of FlatSurf at Point p

Will compute vector even if p is not on surface.
MAKE SURE p IS ON SURFACE

@param p Point to compute normal vector
@param s FlatSurf to compute normal vector of; for s.b > 0 surface posseses translation symmetry along y direction
*/
Vec getNormFlat(Point p, FlatSurf s) {
    double r;
    //if(s.b > 0){
    r = p.x;
    //else{
    //    r = p.y;
    //};
    double den;
    double det = s.k2*s.k2+4*s.k3*(r*r-s.k1);
    if (det > 0){
        den = sqrt(det);
    }
    else{
        return makeVec(-1,0,0); // if the neutron hits the apex of the ellipse we run into a divide by zero problem
    }
    double nx = -2*p.x/den;
    double ny = 0;
    double nz = sign(2*s.k3*p.z+s.k2);
    double n = sqrt(nx*nx+ny*ny+nz*nz);
    return makeVec(nx/n,ny/n,nz/n);
}

/*! \brief Function to compute time of first collision for a ConicSurf

@param p Particle to consider
@param s ConicSurf to consider
@return Time until the propogation or -1 if particle will not hit surface
*/
double getTimeOfFirstCollisionConic(Particle p, ConicSurf s) {
    double tz = (s.zs-p._z)/p._vz;
    if (tz < 0) {
       tz = 0;
       if (p._z > s.ze)
            return -1;
    }

    Particle p2 = copyMoveParticleT(tz,p);

    double A = p2._vx*p2._vx+p2._vy*p2._vy-s.k3*p2._vz*p2._vz;
    double B = 2*(p2._vx*p2._x+p2._vy*p2._y-s.k3*p2._vz*p2._z)-s.k2*p2._vz;
    double C = p2._x*p2._x+p2._y*p2._y-s.k3*p2._z*p2._z-s.k2*p2._z-s.k1;

    double t = solveQuad(A,B,C);

    if (t <= 0 || p2._vz*t+p2._z > s.ze || p2._vz*t+p2._z < s.zs)
        return -1;
    return t+tz;
}

/*! \brief Function to compute time of first collision for a FlatSurf

@param p Particle to consider
@param s FlatSurf to consider
@return Time until the propogation or -1 if particle will not hit surface
*/
//TODO
double getTimeOfFirstCollisionFlat(Particle p, FlatSurf s) {
    double tz = (s.zs-p._z)/p._vz;
    if (tz < 0) {
       tz = 0;
       if (p._z > s.ze)
            return -1;
    }

    Particle p2 = copyMoveParticleT(tz,p);
    double vs = 0;//the vector important for calculating the intersection with the ellipse
    double s0 = 0;
    double vt = 0;//the other component only important for testing whether the mirror is hit
    double t0 = 0;
    //if(s.b > 0){//obsolete iteration allowing to rotate by 90 deg with out rotation in McStas, not really needed
    vs = p2._vx;
    s0 = p2._x;
    vt = p2._vy;
    t0 = p2._y;

    //}
    /*else{
    vs = p2._vy;
    s0 = p2._y;
    vt = p2._vx;
    t0 = p2._x;
    };
    */
    double A = vs*vs-s.k3*p2._vz*p2._vz;
    double B = 2*(vs*s0-s.k3*p2._vz*p2._z)-s.k2*p2._vz;
    double C = s0*s0-s.k3*p2._z*p2._z-s.k2*p2._z-s.k1;

    double t = solveQuad(A,B,C);

    if (t <= 0 || p2._vz*t+p2._z > s.ze || p2._vz*t+p2._z < s.zs||vt*t+t0 < s.ll||vt*t +t0 > s.rl)
        return -1;
    return t+tz;
}

/*! \brief Function to handle supermirror reflectivity copied from mcstas.
@note Uses only m-value for calculating reflectivity curve TODO more sophisticated formulae in the future

@param q k_i - k_f momentum transfer of the neutron at the super mirror surface
@param m supermirror m-value
@param R_0 low angle reflectivity
@param Q_c critical momentum transfer of the super mirror

@return p weight reduction of the neutron for further simulation
*/
double calcSupermirrorReflectivity(double q, double m, double R_0, double Q_c){
    double arg;
    double beta = 0;//values fitting supermirror data from sn
    double alpha = 2.5;
    double W = 0.004;
    double weight = 1.0; //neutron weight to be transformed
    q = fabs(q);
    if (m >= 10){
        weight = 1.0;
        return weight;
    }
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	    alpha=m;
	    beta=0;
        }
    }
    arg = W > 0 ? (q - m*Q_c)/W : 11;
    if (arg > 10 || m <= 0 || Q_c <=0 || R_0 <= 0) {
      weight = 0.0;
      return weight;
    }

    if (m < 1) { Q_c *= m; m=1; }

    if(q <= Q_c) {
      weight = R_0;
      return weight;
    }


    weight = R_0*0.5*(1 - tanh(arg))*(1 - alpha*(q - Q_c) + beta*(q - Q_c)*(q - Q_c));
    return weight;
}

/*! \brief Function to handle reflection of neutron for a ConicSurf.

@note Uses step function for reflectivity

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect
@param s ConicSurf to use

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
double reflectNeutronConic(Particle* p, ConicSurf s) {//TODO add super mirror reflectivity an passing neutrons
    Vec n = getNormConic(getParticlePos(*p),s);
    Vec pv = getParticleVel(*p);
    double weight = 0;
    

    double v = getMagVec(pv);
    double vn = dotVec(pv,n);
    
    weight = calcSupermirrorReflectivity(V2Q_conic*vn*2, s.m, 1.0, 0.0218);

    //Hitting shell from outside

    if (vn > 0 && !s.doubleReflections) {
        absorbParticle(p);
        return -1;
    }



    double ga = fabs(acos(vn/v)) - M_PI/2;
    double gc = 6.84459399932*s.m/v;

    if (weight <= 0) {
        printf("weight <0");
        absorbParticle(p);
        return -1;
    }
    else {
        p->_vx = p->_vx-2*vn*n.x;
        p->_vy = p->_vy-2*vn*n.y;
        p->_vz = p->_vz-2*vn*n.z;
        p->w *= weight;
    }
    return ga;
}

/*! \brief Function to handle refraction of neutron.

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect/refract
@param n Normal vector
@param m1 m-value of the material we come from
@param m2 m-value of the goal material we go to

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
void refractNeutronFlat(Particle* p, Vec n, double m1, double m2) {
    Vec pv = getParticleVel(*p);
    //printf("incoming %.14g %.14g %.14g\n", pv.x, pv.y, pv.z);
    //printf("normal %.14g %.14g %.14g\n", n.x, n.y, n.z);
    double v = getMagVec(pv);
    double vn = dotVec(pv, n);
    //printf("vn %.9g", vn);
    Vec v_p = difVec(pv, skalarVec(n, vn));

    double k2_perp = POT_V*(m1*m1-m2*m2)+vn*vn;
    //printf("the magnitude %f\n", k2_perp);
    if (k2_perp>0){//refraction
        
        k2_perp = sqrt(k2_perp)*sign(vn);
        p->_vx = v_p.x + n.x*k2_perp;
        p->_vy = v_p.y + n.y*k2_perp;
        p->_vz = v_p.z + n.z*k2_perp;
        p-> silicon *= -1;// from silicon to air or vice versa
        //printf("resulting vector %f %f %f\n", p->_vx, p->_vy, p->_vz);
    }else{//total reflection, no change in material
        p->_vx = p->_vx-2*vn*n.x;
        p->_vy = p->_vy-2*vn*n.y;
        p->_vz = p->_vz-2*vn*n.z;
        //no change in material
    }

    
}

/*! \brief Function to handle reflection of neutron for a FlatSurf.

@note Uses step function for reflectivity

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect
@param s FlatSurf to use

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
double reflectNeutronFlat(Particle* p, FlatSurf s) {
    Vec n = getNormFlat(getParticlePos(*p),s);
    Vec pv = getParticleVel(*p);
    //printf("nothing");
    double v = getMagVec(pv);
    double vn = dotVec(pv,n);
    //printf("before %f \n", p->w);
    //Hitting shell from outside For FlatSurface this has to be checked
    // make it able to reflect from the outside
    double ga = fabs(acos(vn/v)) - M_PI/2;
    double gc = 6.84459399932*s.m/v;


    double weight = 0;
    weight = calcSupermirrorReflectivity(V2Q_conic*2*vn, s.m, 0.995, 0.0218);

    if (vn > 0 && !s.doubleReflections) {
        absorbParticle(p);
        return -1;
    }


    if (weight < 0) {
        printf("this happens?");
        absorbParticle(p);
        return -1;
    } else {//here we need to implement the refraction
        if (getRandom() <= weight){//to be updated to use the mcstas random function or quasoi deterministic model
            //printf("oh a reflections\n");
            //printf("this total reflection?");
            p->_vx = p->_vx-2*vn*n.x;
            p->_vy = p->_vy-2*vn*n.y;
            p->_vz = p->_vz-2*vn*n.z;
            return ga;
        }
        else{//
            //printf("enter the refraction process");
            double m1 = 0;
            double m2 = 0;
            //if no mirrorwidth is specified no refraction has to be calc
            if (p->silicon == 0){
                return ga;
            }
            if (p->silicon == 1){
                m1 = m_Si;
                m2 = 0;
            }
            if (p->silicon == -1){
                m1 = 0;
                m2 = m_Si;                
            }
            //printf("k1 %f k2 %f k3 %f x %f y %f z %f\n", s.k1, s.k2, s.k3, p->_x, p->_y, p->_z);
            refractNeutronFlat(p, n, m1, m2);// this can still lead to total reflection, we miss the supermirror, but are still reflected by the silicon takes care of change of material for refraction
            return ga;
        }
    }
    return ga;
}


/*! \brief Function to handle raytracing of neutron for a ConicSurf.

@param p Pointer of particle to reflect
@param c ConicSurf to use

*/
void traceNeutronConic(Particle* p, ConicSurf c) {
    double t = getTimeOfFirstCollisionConic(*p, c);
    if (t < 0)
        return;
    else {
        moveParticleT(t, p);
        double ga = reflectNeutronConic(p, c);
#if REC_MAX_GA
        if (ga > c.max_ga) {
            c.max_ga = ga;
            c.max_ga_z0 = p->_z;
        }
#endif
    }
}

/*! \brief Function to handle raytracing of neutron for a FlatSurf.

@param p Pointer of particle to reflect
@param f FlatSurf to use

*/
void traceNeutronFlat(Particle* p, FlatSurf f) {
    double t = getTimeOfFirstCollisionFlat(*p, f);
    if (t < 0)
        return;
    else {

        moveParticleT(t, p);

        //printf("weight before reflect %f", p->w);
        double ga = reflectNeutronFlat(p, f);
        //printf("weight after reflect %f\n", p->w);
#if REC_MAX_GA
        if (ga > f.max_ga) {
            f.max_ga = ga;
            f.max_ga_z0 = p->_z;
        }
#endif
    }
}
/** @} */ //end of conicgroup
/////////////////////////////////////
// Scene Functions
/////////////////////////////////////
/** @ingroup simgroup
    @{
*/
enum GEO {
    NONE,
    DETECTOR,
    DISK,
    CONIC,
    FLAT
};

//! Function to generate an empty Scene
Scene makeScene() {
    Scene s;
    s.num_f = 0;
    s.num_c = 0;
    s.num_di = 0;
    s.num_d = 0;

    s.traceNeutronFlat = traceNeutronFlat;
    s.traceNeutronConic = traceNeutronConic;
    s.traceNeutronDisk = traceNeutronDisk;
    s.traceNeutronDetector = traceNeutronDetector;

    return s;
}

//! Function to init simulation items
/*! Should be called after all items
have been added to scene but before
neutrons are traced.

@param s Pointer of Scene to init
*/
void initSimulation(Scene* s) {
    //
}

/*! \brief Function to raytrace single neutron through geometries specified by d, di and c.

@param p Pointer of particle to trace
@param s Scene to trace
*/
void traceSingleNeutron(Particle* p, Scene s) {

    int contact = 1;
    do {
        double t;
        enum  GEO type = NONE;
        int index = -1;
        int i;

        for (i = 0; i < s.num_c; i++) {
            double t2 = getTimeOfFirstCollisionConic(*p,s.c[i]);

            if (t2 <= 0)
                continue;
            if (index == -1 || t2 < t) {
                type = CONIC;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_f; i++) {
            double t2 = getTimeOfFirstCollisionFlat(*p,s.f[i]);

            if (t2 <= 0)
                continue;
            if (index == -1 || t2 < t) {
                type = FLAT;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_di; i++)  {
            double t2 = getTimeOfFirstCollisionDisk(*p,s.di[i]);

            if (t2 <= 0)
                continue;
            else if (index == -1 || t2 < t) {
                type = DISK;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_d; i++) {
            double t2 = getTimeOfFirstCollisionDetector(*p,s.d[i]);

            if (t2 <= 0)
                continue;
            else if (index == -1 || t2 < t) {
                type = DETECTOR;
                index = i;
                t = t2;
            }
        }

        switch (type) {
            case DETECTOR:
                s.traceNeutronDetector(p, s.d[index]);
                break;
            case FLAT:
                s.traceNeutronFlat(p, s.f[index]);
                break;
            case DISK:
                s.traceNeutronDisk(p, s.di[index]);
                break;
            case CONIC:
                s.traceNeutronConic(p, s.c[index]);
                break;
            default:
                contact = 0;
                break;
        }
    } while (contact && !p->absorb);

}

//!Finishes tracing the scene
/*! This function should be called after all of the
particles have been raytraced.

@param s Pointer of Scene to finish tracing
*/
void finishSimulation(Scene* s) {
    int i;

    //Finish Detectors
    for (i=0; i < s->num_d; i++)
        finishDetector(s->d[i]);
}

/** @} */ //end of ingroup simgroup

#endif
