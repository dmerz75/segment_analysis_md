// readfile.h
#ifndef _CURVATURE_H_
#define _CURVATURE_H_

// libraries:
/* #include <stdio.h> */
/* #include <stdlib.h> */
/* #include <string> */
#include <vector>

// Definitions:
/* #define BUFFERSIZE 900 */


/* ---------------------------------------------------------
   Class declarations
   --------------------------------------------------------- */
/* #ifndef CHAIN_H */
/* #define CHAIN_H */
#include "chain.h"
/* #endif */

class Curvature {

public:
    Vector p1;
    Vector p2;
    Vector p3;

    double radius;
    double curvature;

    Curvature(); // Constructor declared.

    // functions:
    void get_curvature();
    void print_points();
    void print_rc();
    void points_ascending_z();
    void get_circumscribed_circle();
    void get_radius_circumscribed();
};
inline Curvature::Curvature() {
    // constructor
    radius = 0.0;
    curvature = 0.0;
    /* p1.x = 0.0; */
    /* p1.y = 0.0; */
    /* p1.z = 0.0; */
}
inline void Curvature::print_points() {
    printf("p1: %.6f %.6f %.6f\n",p1.x,p1.y,p1.z);
    printf("p2: %.6f %.6f %.6f\n",p2.x,p2.y,p2.z);
    printf("p3: %.6f %.6f %.6f\n",p3.x,p3.y,p3.z);
}
inline void Curvature::print_rc() {
    printf("the radius: %f\n",radius);
    printf("the curvature: %f\n",curvature);
}
inline void Curvature::get_circumscribed_circle() {

    /* if  */
    /* http://stackoverflow.com/questions/13977354/build-circle-from-3-points-in-3d-space-implementation-in-c-or-c */

    // triangle "edges"
    /* const Vector3d t = p2-p1; */
    /* const Vector3d u = p3-p1; */
    /* const Vector3d v = p3-p2; */


    Vector t, u, v, w, circCenter, circAxis;
    double circRadius;

    circCenter.x = circCenter.y = circCenter.z = 0.0;
    circAxis.x = circAxis.y = circAxis.z = 0.0;
    circRadius = 0.0;

    radius = circRadius;

    if ((isnan(p1.x) == 1) || (isnan(p1.y) == 1) || (isnan(p1.z) == 1))
    {
        debug("p1 is nan.\n");
        return;
    }
    if ((isnan(p2.x) == 1) || (isnan(p2.y) == 1) || (isnan(p2.z) == 1))
    {
        debug("p2 is nan.\n");
        return;
    }
    if ((isnan(p3.x) == 1) || (isnan(p3.y) == 1) || (isnan(p3.z) == 1))
    {
        debug("p3 is nan.\n");
        return;
    }


    t = get_vector(p1,p2);
    u = get_vector(p1,p3);
    v = get_vector(p2,p3);
    /* w = cross_product(t,v); */
    w = cross_product(t,u);

#ifndef NDEBUG
    t.print_Vector();
    u.print_Vector();
    v.print_Vector();
    w.print_Vector();
#endif

    const double wsl = magnitude(w);

    /* // triangle normal */
    /* const Vector3d w = t.crossProduct(u); */
    /* const double wsl = w.getSqrLength(); */

    /* // area of the triangle is too small (you may additionally check the points */
    /* // for colinearity if you are paranoid) */

    /* if (wsl<10e-14) return false; */
    /* if (wsl<10e-14) exit(1); */

    /* // helpers */
    const double iwsl2 = 1.0 / (2.0*wsl);
    /* const double tt = t*t; */
    /* const double uu = u*u; */


    if (wsl<10e-14)
    {
        return;
    }
    if (isnan(wsl))
    {
        debug("wsl is nan!\n");
        return;
    }
    if (isnan(iwsl2))
    {
        debug("wsl is a num, iwsl2 is nan!\n");
        return;
    }


    debug("wsl: %f   iwsl2: %f\n",wsl,iwsl2);

    double tt, uu, vv, oneover_sqrtwsl;
    tt = dot_product(t,t);
    uu = dot_product(u,u);
    vv = dot_product(v,v);


    /* // result circle */
    /* Vector3d circCenter = p1 + (u*tt*(u*v) - t*uu*(t*v)) * iwsl2; */
    /* double   circRadius = sqrt(tt * uu * (v*v) * iwsl2*0.5); */

    /* Vector uv,tv,utt,tuu,uttuv,tuutv,dif,prod; */
    /* Vector utt,tuu,uttuv,tuutv,dif,prod; */
    Vector uttuv,tuutv,dif,prod;
    double uv,tv,ttuv,uutv;
    uv = dot_product(u,v);
    tv = dot_product(t,v);
    /* utt = scalar_mult(u,tt); */
    ttuv = tt * uv;
    uutv = uu * tv;
    /* tuu = scalar_mult(t,uu); */
    /* uttuv = dot_product(utt,uv); */
    /* tuutv = dot_product(tuu,tv); */
    uttuv = scalar_mult(u,ttuv);
    tuutv = scalar_mult(t,uutv);

    dif = difference(tuutv,uttuv);
    prod = scalar_mult(dif,iwsl2);
    circCenter = vec_add(p1,prod);
    circRadius = sqrt(tt * uu * vv * iwsl2 * 0.5);


    /* Vector3d circAxis   = w / sqrt(wsl); */
    oneover_sqrtwsl = 1 / sqrt(wsl);
    circAxis = scalar_mult(w,oneover_sqrtwsl);

    debug("the Center/Axis | Radius are: %.3f %.3f %.3f, %.3f %.3f %.3f  |  %f\n", \
           circCenter.x,circCenter.y,circCenter.z, \
           circAxis.x,circAxis.y,circAxis.z, \
           circRadius);

    radius = circRadius;
    curvature = abs( 1 / circRadius);

}
inline void Curvature::get_radius_circumscribed()
{
    /* https://www.physicsforums.com/threads/equation-of-a-circle-through-3-points-in-3d-space.173847/ */

    Vector a, b, c;
    a = get_vector(p1,p2);
    b = get_vector(p2,p3);
    c = get_vector(p3,p1);

    // magnitude
    double am, bm, cm;

    am = magnitude(a);
    bm = magnitude(b);
    cm = magnitude(c);

    double abc, ab2, bc2, ca2, a4, b4, c4, under, denom;
    abc = ab2 = bc2 = ca2 = a4 = b4 = c4 = under = denom = 0.0;

    abc = am * bm * cm;
    ab2 = 2 * am * am * bm * bm;
    bc2 = 2 * bm * bm * cm * cm;
    ca2 = 2 * cm * cm * am * am;
    a4  = am * am * am * am * -1;
    b4  = bm * bm * bm * bm * -1;
    c4  = cm * cm * cm * cm * -1;

    under = ab2 + bc2 + ca2 + a4 + b4 + c4;

    if (under > 0.0)
    {
        denom = sqrt(under);
        radius = abc / denom;
    }
    else
    {
        debug("negative sqrt(denom)\n");
    }
}
inline void Curvature::points_ascending_z() {
    /* print_points(); */

    Vector x1,x2,x3;
    x1.x = p1.x;
    x1.y = p1.y;
    x1.z = p1.z;
    x2.x = p2.x;
    x2.y = p2.y;
    x2.z = p2.z;
    x3.x = p3.x;
    x3.y = p3.y;
    x3.z = p3.z;


    // find lowest.
    if (x1.z < x2.z)
    {
        if (x1.z < x3.z) // x1.z is lowest: 1 < 2, 1 < 3
        {
            // p1 stays.
            if (x2.z < x3.z)
            {
                // nothing. 1 < 2 < 3.
            }
            else
            {
                // 1 < 2, 1 < 3, 3 < 2; switch 2 to 3, 3 to 2.
                p2.x = x3.x;
                p2.y = x3.y;
                p2.z = x3.z;
                p3.x = x2.x;
                p3.y = x2.y;
                p3.z = x2.z;
            }
        }
        else
        {
            // 1 < 2, 3 < 1; 3 < 1 < 2.
            // p1-> x3  p2-> x1  p3-> x2
            p1.x = x3.x;
            p1.y = x3.y;
            p1.z = x3.z;
            p2.x = x1.x;
            p2.y = x1.y;
            p2.z = x1.z;
            p3.x = x2.x;
            p3.y = x2.y;
            p3.z = x2.z;
        }
    }
    else
    {
        // 2 < 1
        if (x1.z < x3.z)
        {
            // 2 < 1, 1 < 3;  2 < 1 < 3;
            // p1->2, p2->1, p3->3
            p1.x = x2.x;
            p1.y = x2.y;
            p1.z = x2.z;
            p2.x = x1.x;
            p2.y = x1.y;
            p2.z = x1.z;
        }
        else
        {
            // 2 < 1, 3 < 1; 3?2 < 1
            if (x2.z < x3.z)
            {
                // 2 < 1, 3 < 1, 2 < 3; 3?2 < 1, + 2 < 3 ... 2 < 3 < 1;
                // p1->2, p2->3, p3->1
                p1.x = x2.x;
                p1.y = x2.y;
                p1.z = x2.z;
                p2.x = x3.x;
                p2.y = x3.y;
                p2.z = x3.z;
                p3.x = x1.x;
                p3.y = x1.y;
                p3.z = x1.z;
            }
            else
            {
                // 2 < 1, 3 < 1, 3?2 < 1, 3 < 2 ----> 3 < 2 < 1
                // p1->3, p2->2, p3->1
                p1.x = x3.x;
                p1.y = x3.y;
                p1.z = x3.z;
                p3.x = x1.x;
                p3.y = x1.y;
                p3.z = x1.z;
            }
        }
    }


    /* // sorting! */
    /* if ((x1.z <= x2.z) && (x1.z <= x3.z)) // 1 is already least. */
    /* { */
    /*     if (x3.z < x2.z) // 1 is least, 3 < 2, switch 2 & 3. */
    /*     { */
    /*         p2.x = x3.x; */
    /*         p2.y = x3.y; */
    /*         p2.z = x3.z; */
    /*         p3.x = x2.x; */
    /*         p3.y = x2.y; */
    /*         p3.z = x2.z; */
    /*     } */
    /* } */
    /* else if (x2.z <= x3.z) // 1 is not least, 2 is least. p1 = 2. */
    /* { */
    /*     p1.x = x2.x; */
    /*     p1.y = x2.y; */
    /*     p1.z = x2.z; */

    /*     if (x1.z <= x3.z) // 1 < 3. p2 = 1, p3 stays 3 */
    /*     { */
    /*         p2.x = x1.x; */
    /*         p2.y = x1.y; */
    /*         p2.z = x1.z; */
    /*     } */
    /*     else // 1 > 3. p2 = 3, p3 = 1. */
    /*     { */
    /*         // 3 is 2nd, 1 is greatest. */
    /*         p2.x = x3.x; */
    /*         p2.y = x3.y; */
    /*         p2.z = x3.z; */
    /*         p3.x = x1.x; */
    /*         p3.y = x1.y; */
    /*         p3.z = x1.z; */
    /*     } */
    /* } */
    /* else // 3 is least. */
    /* { */
    /*     // p3 stays x3 */
    /*     p1.x = x3.x; */
    /*     p1.y = x3.y; */
    /*     p1.z = x3.z; */

    /*     if (x1.z > x2.z) // 2 is least. 1 is mid. */
    /*     { */
    /*         p2.x = x1.x; */
    /*         p2.y = x1.y; */
    /*         p2.z = x1.z; */
    /*         p1.x = x2.x; */
    /*         p1.y = x2.y; */
    /*         p1.z = x2.z; */
    /*     } */
    /*     else // 3-least, 1 is 2nd, 2 is greatest */
    /*     { */
    /*         p3.x = x2.x; */
    /*         p3.y = x2.y; */
    /*         p3.z = x2.z; */
    /*         p2.x = x1.x; */
    /*         p2.y = x1.y; */
    /*         p2.z = x1.z; */
    /*     } */
    /* } */

    if ((p1.z > p2.z) || (p1.z > p3.z))
    {
        printf("points_ascending_z has failed. p1.z is greater than p2.z or p3.z\n");
        exit(1);
    }
    if (p2.z > p3.z)
    {
        printf("points_ascending_z has failed. p2.z > p3.z\n");
        exit(1);
    }
}
inline void Curvature::get_curvature() {
    // from 3 points
    // http://www.intmath.com/applications-differentiation/8-radius-curvature.php
    // z = 0.

    double m1,m2,xc,yc;
    // get_vector(p1,p2,&m1);
    // get_vector(p1,p3,&m2);
    m1 = (p2.y - p1.y) / (p2.x - p1.x);
    m2 = (p3.y - p2.y) / (p3.x - p2.x);

    double mmyy, mxx2, mxx3, mm2;
    mmyy = m1 * m2 * ( p1.y - p3.y );
    mxx2 = m2 * (p1.x + p2.x);
    mxx3 = -1.0 * m1 * (p2.x + p3.x);
    mm2  = 2 * (m2 - m1);
    xc = ( mmyy + mxx2 + mxx3 ) / mm2;

    double minus_onem, xx2, yy2;
    minus_onem = -1.0 / m1;
    xx2 = (p1.x + p2.x) / 2;
    yy2 = (p1.y + p2.y) / 2;
    yc = minus_onem * (xc - xx2) + yy2;

    radius = sqrt(pow((p2.x - xc),2) + pow((p2.y - yc),2));

    // p1-p3 / 2  is opposite.
    // radius is hypotenuse.
    double d1,sintheta;
    d1 = sintheta = 0.0;
    d1 = distance(p1,p3);

    sintheta = d1/( 2 * radius);
    curvature = 180.0 * 2 * asin(sintheta) / M_PI;




    /* if (radius > 2000.0) { */
    /*     radius = 2000.01; */
    /* } */

    /* // isnan returns 1 for true, as in nan. */
    /* if (isnan(curvature) != 0 ) { */
    /*     curvature = 0.0; */
    /* } else if (curvature > 500.0) { */
    /*     curvature = 500.01; */
    /* } */

    debug("xc: %f\n",xc);
    debug("yc: %f\n",yc);
    debug("radius: %f\n",radius);
    debug("curvature: %f\n",curvature);
}



/* ---------------------------------------------------------
   function declarations
   --------------------------------------------------------- */
/* Curvature compute_curvature(Chain *chain,int i,FILE *fp_curv, FILE *fp_rad); */
Curvature compute_curvature3(Chain *chain,int p1,int p2,int p3,FILE *fp_curv, FILE *fp_rad);
Curvature compute_curvature4(Chain *chain,int p1,int p2,int p3,int p4,FILE *fp_curv, FILE *fp_rad);
Vector translate_to_origin(Chain *chain,int num_chain,int max_num);
Matrix rotate_system_around_axis(Chain *chain,int max_num_chains,       \
                                 int axisfrom,int axisto,Vector Origin,Vector Endpoint,int Dir);

// testing:
/* std::vector<int> get_protofilaments(Chain *chain,int chains_to_use); */
std::vector<int> get_monomers_with_no_southern_neighbor(Chain *chain,int chains_to_use);
std::vector<int> get_protofilament(Chain *chain,int pf);

// probably not that useful ..
void get_xy_plane(Chain *chain,int max_num_chains,int num_chain);
void rotate_system_around_z_to_y(Chain *chain,int begin,int end);
Vector midpoint_of_indices_of_2_chains(Chain *chain,int c1,int c2,int &idx1, int &idx2);
int get_pos_closest_to_point(Chain *chain,int c,Vector v);

#endif
