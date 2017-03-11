// md.h
#ifndef _MD_H_
#define _MD_H_

#include <stdio.h>
#include <math.h>
#include "debug.h"


/* !!!!!!!!!!! */
/* struct Vector { */
/*     double x; */
/*     double y; */
/*     double z; */
/* }; */

/* class dvalue { */
/* public: */
/*     double d[32]; */
/*     dvalue (); // Constructor declared. */
/* }; */

/* inline dvalue::dvalue() { */
/*     for (int i=0; i<32; i++) { */
/*         d[i] = 0.0; */
/*     } */
/* } */


class Vector {
public:
    double x;
    double y;
    double z;
    Vector (); // Constructor declared.
    void print_Vector();
};
inline Vector::Vector() {
    x = 0.0;
    y = 0.0;
    z = 0.0;
}
inline void Vector::print_Vector() {
    printf("\nxyz: %f %f %f\n",x,y,z);
}

class Matrix {
public:
    double x1,x2,x3;
    double y1,y2,y3;
    double z1,z2,z3;

    // Constructor
    Matrix (); // Constructor declared.

    // Destructor
    ~Matrix() {};

    // declarations
    /* void build_rotation(double theta,char sel[1]); */
    /* void build_rotation(double theta,string sel); */
    void build_rotation(double theta,int sel);
    void transpose();
    void print_Matrix();


    /* double theta; */
    /* char sel; */
};
inline Matrix::Matrix() {
    x1 = x2 = x3 = 0.0;
    y1 = y2 = y3 = 0.0;
    z1 = z2 = z3 = 0.0;
}
inline void Matrix::build_rotation(double theta,int sel) {
    // x:0, y:1, z:2;
    debug("theta: %f\n",theta);
    debug("selection: %d\n",sel);

    double costheta,sintheta,msintheta;
    costheta = cos(theta);
    sintheta = sin(theta);
    msintheta = -1.0 * sintheta;
    if(sel == 0) { // x
        x1 = 1.0;
        x2 = x3 = 0.0;
        y1 = z1 = 0.0;
        y2 = z3 = costheta;
        y3 = msintheta;
        z2 = sintheta;
    } else if (sel == 1) { // y
        y2 = 1.0;
        y1 = y3 = 0.0;
        x2 = z2 = 0.0;
        x1 = z3 = costheta;
        x3 = msintheta;
        z1 = sintheta;
    } else if (sel == 2) { // z
        z3 = 1.0;
        z1 = z2 = 0.0;
        x3 = y3 = 0.0;
        x1 = y2 = costheta;
        x2 = msintheta;
        y1 = sintheta;
    }
}
inline void Matrix::print_Matrix() {
    // x:0, y:1, z:2;
    printf("X: %f %f %f\n",x1,x2,x3);
    printf("Y: %f %f %f\n",y1,y2,y3);
    printf("Z: %f %f %f\n",z1,z2,z3);
}
inline void Matrix::transpose() {
    // x:0, y:1, z:2;
    // x1,y2,z3;
    double tempx2,tempx3,tempy3;
    tempx2 = x2;
    tempx3 = x3;
    tempy3 = y3;
    // upper
    x2 = y1;
    x3 = z1;
    y3 = z2;
    // lower
    y1 = tempx2;
    z1 = tempx3;
    z2 = tempy3;
}


// declarations
double distance ( Vector v1, Vector v2 );
Vector difference ( Vector v1, Vector v2 );
double magnitude ( Vector v1 );
double get_costheta ( Vector v1, Vector v2 );
double dot_product ( Vector v1, Vector v2 );
double get_sintheta ( Vector v1, Vector v2 );
/* Vector midpoint ( Vector v1 ); */
Vector midpoint2(Vector v1,Vector v2);
Vector cross_product ( Vector v1, Vector v2 );
Vector scalar_mult ( Vector v1, double a );
Vector vec_add ( Vector v1, Vector v2 );
Vector prod_matrix_vector(Matrix M, Vector V);
/* void print_vector(Vector V); */

/* double curvature(Vector p1, Vector p2, Vector p3); */

/* can't do it by return type, Yes for variable input */
// overloaded declarations
Vector get_vector (Vector v1,Vector v2,Vector *v3 );
Vector get_vector (Vector v1,Vector v2);
Vector normalize (Vector v1,Vector *v2);
Vector normalize (Vector v1);

#endif
