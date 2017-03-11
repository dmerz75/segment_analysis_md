// md.c

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
#include <math.h>

/* ---------------------------------------------------------
   other headers
   --------------------------------------------------------- */
#include "md.h"



/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
double distance ( Vector v1, Vector v2 ) {
    /* Distance between 2 points! */
    double distance;
    double xdist, ydist, zdist;

    xdist = v2.x - v1.x;
    ydist = v2.y - v1.y;
    zdist = v2.z - v1.z;

    distance = sqrt( xdist * xdist + ydist * ydist + zdist * zdist);
    return distance;
}
Vector difference ( Vector v1, Vector v2 ) {
    /* difference */
    Vector v;

    v.x = v2.x - v1.x;
    v.y = v2.y - v1.y;
    v.z = v2.z - v1.z;

    return v;
}
double magnitude ( Vector v1 ) {
    /* get magnitude */
    double magnitude;
    magnitude = sqrt( v1.x * v1.x + v1.y * v1.y + v1.z * v1.z );
    return magnitude;
}
// Vector midpoint ( Vector v1 ) {
//     /* get midpoint */
//     Vector v2;
//     // magnitude = sqrt( v1.x * v1.x + v1.y * v1.y + v1.z * v1.z );
//     v2.x = v1.x * 0.5;
//     v2.y = v1.y * 0.5;
//     v2.z = v1.z * 0.5;
//     return v2;
// }
Vector midpoint2(Vector v1,Vector v2){
    /* get midpoint */
    Vector v;
    v.x = (v1.x + v2.x) * 0.5;
    v.y = (v1.y + v2.y) * 0.5;
    v.z = (v1.z + v2.z) * 0.5;
    return v;
}
double get_costheta ( Vector v1, Vector v2 ) {
    /* !!! v1 & v2 need to be normalized. */
    double costheta;
    costheta = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    return costheta;
}
double dot_product ( Vector v1, Vector v2 ) {
    /* !!! v1 & v2 not normalized, then it's a dot product. */
    double dotp;
    dotp = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    return dotp;
}
double get_sintheta ( Vector v1, Vector v2 ) {
    /* !!! v1 & v2 need to be normalized. */
    Vector v;
    double sintheta;
    // double m;
    // costheta = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = (v1.x * v2.z - v1.z * v2.x ) * -1;
    v.z = v1.x * v2.y - v1.y * v2.x;
    sintheta = magnitude(v);
    return sintheta;
}
Vector cross_product ( Vector v1, Vector v2 ) {
    // cross product
    // row1: i j k
    // row2: v1
    // row3: v2
    Vector v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = (v1.x * v2.z - v1.z * v2.x ) * -1;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}
Vector scalar_mult ( Vector v1, double a ){
    // scalar_mult
    Vector v;
    v.x = v1.x * a;
    v.y = v1.y * a;
    v.z = v1.z * a;
    return v;
}
Vector vec_add ( Vector v1, Vector v2 ){
    // scalar_mult
    Vector v;
    v.x = v1.x + v2.x;
    v.y = v1.y + v2.y;
    v.z = v1.z + v1.z;
    return v;
}
Vector prod_matrix_vector(Matrix M, Vector V) {
    Vector v;
    v.x = M.x1 * V.x + M.x2 * V.y + M.x3 * V.z;
    v.y = M.y1 * V.x + M.y2 * V.y + M.y3 * V.z;
    v.z = M.z1 * V.x + M.z2 * V.y + M.z3 * V.z;

    return v;
}
// void print_Vector(Vector V){
//     printf("print_Vector(xyz): %f %f %f\n",V.x,V.y,V.z);
// }

/* ---------------------------------------------------------
   overloaded functions
   --------------------------------------------------------- */
Vector normalize (Vector v1, Vector *v2) {
    /* v1 is vector, v2 is normalized, destination */
    double magnitude;
    Vector v2b;
    magnitude = sqrt( v1.x * v1.x + v1.y * v1.y + v1.z * v1.z );

    v2b.x = v2->x = v1.x / magnitude;
    v2b.y = v2->y = v1.y / magnitude;
    v2b.z = v2->z = v1.z / magnitude;

    return v2b;
}
Vector normalize (Vector v1) {
    /* v1 is vector, v2 is normalized, destination */
    double magnitude;
    Vector v2b;
    magnitude = sqrt( v1.x * v1.x + v1.y * v1.y + v1.z * v1.z );

    v2b.x = v1.x / magnitude;
    v2b.y = v1.y / magnitude;
    v2b.z = v1.z / magnitude;

    return v2b;
}
Vector get_vector ( Vector v1, Vector v2, Vector *v3 ) {
    /* vector between 2 points! */
    Vector v;
    v.x = v3->x = v2.x - v1.x;
    v.y = v3->y = v2.y - v1.y;
    v.z = v3->z = v2.z - v1.z;
    return v;
}
Vector get_vector ( Vector v1, Vector v2) {
    /* vector between 2 points! */
    Vector v;
    v.x = v2.x - v1.x;
    v.y = v2.y - v1.y;
    v.z = v2.z - v1.z;
    return v;
}
