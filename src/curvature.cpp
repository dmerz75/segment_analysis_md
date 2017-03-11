// curvature.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h> // strtod?, stod
#include <math.h> // M_PI
#include <vector>

// namespaces
// using namespace

// extern
extern int endoffile;

/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "readfile.h"
#include "chain.h"
#include "md.h"
#include "curvature.h"


/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
// Curvature compute_curvature(Chain *chain,int i,FILE *fp_curv,FILE *fp_rad) {
Curvature compute_curvature3(Chain *chain,int p1,int p2, int p3,FILE *fp_curv,FILE *fp_rad) {
    // curvature was moved to //md.cpp// curvature.h
    Chain *chain_local;
    try {
        chain_local = new Chain [3];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        exit(1);
        // return 1;
    }
    // int p1,p2,p3,p4;

    // COPY num_atoms_ca, and the pos x,y,z
    chain_local[0].num_atoms_ca = chain[p1].num_atoms_ca;
    memcpy(chain_local[0].pos,chain[p1].pos,sizeof(chain[p1].pos));
    chain_local[1].num_atoms_ca = chain[p2].num_atoms_ca;
    memcpy(chain_local[1].pos,chain[p2].pos,sizeof(chain[p2].pos));
    chain_local[2].num_atoms_ca = chain[p3].num_atoms_ca;
    memcpy(chain_local[2].pos,chain[p3].pos,sizeof(chain[p3].pos));

    debug("Computing local centroid.\n");
    chain_local[0].ComputeCentroid();
    chain_local[1].ComputeCentroid();
    chain_local[2].ComputeCentroid();


#ifndef NDEBUG
    debug("original:\n");
    chain[0].print_centroid();
    chain[1].print_centroid();
    chain[2].print_centroid();
    // Compute the Centroids.
    debug("local:\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    // exit(0);
#endif

    // Rotation. (locally)
    Vector t0;
    // ,t1;
    Matrix R1,R2,R3;
    // Matrix R21,R22,R23;
    // R4;
    // Matrix yx01,x12;
    double spot;

    Curvature curv3;
    // Curvature curvright;
    // Curvature curvavg;

    // translate to origin
    t0 = translate_to_origin(chain_local,0,3);
    // printf("after local curvature translation:\n");
    // printf("local:\n");
    // chain_local[0].print_centroid();
    // chain_local[1].print_centroid();
    // chain_local[2].print_centroid();
    // chain_local[3].print_centroid();
    // printf("original:\n");
    // chain[0].print_centroid();
    // chain[1].print_centroid();
    // chain[2].print_centroid();
    // chain[3].print_centroid();
    // exit(0);


    // around Z to X axis, more X.
    // axisfrom Z(2) (no Z comp.), axisto X(0) (UNIT.X); sent to X.
    spot = chain_local[1].centroid.y;
    R1 = rotate_system_around_axis(chain_local,3,2,0,           \
                                   chain_local[0].centroid,     \
                                   chain_local[1].centroid,1);

#ifndef NDEBUG
    debug("pre-round 1.\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
#endif

    // printf("y-spot vs. new y:   %f    %f\n",spot,chain_local[1].centroid.y);
    // Check1.
    if (fabs(chain_local[1].centroid.y) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 1 !!!\n");

        // TRANSPOSE
        R1.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<3; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R1,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

        // rotate twice.
        for(int i=0; i<3; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R1,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();

        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        debug("WARNING: had to exit at compute_curvature/1st rotation.\n");
#endif
        // exit(1);
    } // end of Check1.
    // exit(0);
    debug("----------Rotate #1 Done.-----------\n");



    spot = chain_local[1].centroid.z;
    R2 = rotate_system_around_axis(chain_local,3,1,0,           \
                                   chain_local[0].centroid,     \
                                   chain_local[1].centroid,-1);

#ifndef NDEBUG
    debug("pre-round 2.\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
#endif

    // printf("y-spot vs. new y:   %f    %f\n",spot,chain_local[1].centroid.y);
    // Check1.
    if (fabs(chain_local[1].centroid.z) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 2 !!!\n");

        // TRANSPOSE
        R2.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<3; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R2,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R2,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();


        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();

        debug("WARNING: had to exit at compute_curvature/2nd rotation.\n");
#endif
    } // end of Check1.
    debug("----------Rotate #2 Done.-----------\n");
    // exit(0);


    spot = chain_local[2].centroid.z;
    R3 = rotate_system_around_axis(chain_local,3,0,1,          \
                                    chain_local[1].centroid,    \
                                    chain_local[2].centroid,1);


#ifndef NDEBUG
    debug("pre-round 3.\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();

#endif


    // printf("z-spot vs. new z:   %f    %f\n",spot,chain_local[2].centroid.z);
    // Check1.
    if (fabs(chain_local[2].centroid.z) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 3 !!!\n");

        // TRANSPOSE
        R3.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R3,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R3,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();
        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        debug("WARNING: had to exit at compute_curvature/3rd rotation.\n");
#endif
    } // end of Check1.
    debug("----------Rotate #3 Done.-----------\n");


#ifndef NDEBUG
    // print centroid.
    debug("curvature_check:\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
#endif
    // quality check.
    printf("after_first_rotation(pos-2):\n");
    chain_local[2].print_centroid();

    curv3.p1 = chain_local[0].centroid;
    curv3.p2 = chain_local[1].centroid;
    curv3.p3 = chain_local[2].centroid;
    curv3.get_curvature();

#ifndef NDEBUG
    printf("curv left only, before:\n");
    curv3.print_rc();
#endif

//     // get average
//     if ((curv3.radius < 2000.0) && (curvright.radius < 2000.0)) {
//         curvavg.radius = (curv3.radius + curvright.radius) * 0.5;
//     } else if ((curv3.radius > 2000.0) || (curvright.radius > 2000.0)) {
//         curvavg.radius = 2000.01;
//     }

//     if ((curv3.curvature < 500.0) && (curvright.curvature < 500.0)) {
//         curvavg.curvature = (curv3.curvature + curvright.curvature) * 0.5;
//     } else if ((curv3.curvature > 500.0) || curvright.curvature > 500.0) {
//         curvavg.curvature = 500.01;
//     }


// #ifndef NDEBUG
//     debug("average:\n");
//     curvavg.print_rc();
// #endif

// #ifndef NDEBUG
//     debug("difference:\n");
//     debug("\tradius: %f \n\tcurvature: %f\n",curv3.radius - curvavg.radius,curv3.curvature - curvavg.curvature);
// #endif

    fprintf(fp_curv,"%9.4f ",curv3.curvature);
    fprintf(fp_curv,"%9.4f ",curv3.p3.z);

    fprintf(fp_rad,"%9.4f ",curv3.radius);
    fprintf(fp_rad,"%9.4f ",curv3.p3.z);

    // delete
    delete [] chain_local;
    return curv3;
} // End of compute_curvature3

Curvature compute_curvature4(Chain *chain,int p1,int p2, int p3, int p4,FILE *fp_curv,FILE *fp_rad) {
    // curvature was moved to //md.cpp// curvature.h

    // rotate locally??
    Chain *chain_local;
    try {
        chain_local = new Chain [4];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        exit(1);
        // return 1;
    }

    // int p1,p2,p3,p4;
    // // MT: LonS,i,LonN,LonN-LonN
    // // PF:
    // p1 = chain[i].LonS;
    // p2 = i;
    // p3 = chain[i].LonN;

    // debug("neighbors: %d %d %d\n",p1,p2,p3);
    // if ((p1 == -1) && (p3 == -1)) {
    //     // printf("switch to PF\n");
    //     p1 = i - 1;
    //     p3 = i + 1;
    //     p4 = i + 2;
    //     debug("PF neighbors: %d %d %d %d\n",p1,p2,p3,p4);
    // } else {
    //     debug("MT neighbors: %d %d %d %d\n",p1,p2,p3,p4);

    //     // COULD BE A MAJOR CAUSE OF SEGFAULT
    //     // if chain[i].LonN does not exist.
    //     p4 = chain[chain[i].LonN].LonN;
    // }


    // COPY num_atoms_ca, and the pos x,y,z
    chain_local[0].num_atoms_ca = chain[p1].num_atoms_ca;
    memcpy(chain_local[0].pos,chain[p1].pos,sizeof(chain[p1].pos));
    chain_local[1].num_atoms_ca = chain[p2].num_atoms_ca;
    memcpy(chain_local[1].pos,chain[p2].pos,sizeof(chain[p2].pos));
    chain_local[2].num_atoms_ca = chain[p3].num_atoms_ca;
    memcpy(chain_local[2].pos,chain[p3].pos,sizeof(chain[p3].pos));
    chain_local[3].num_atoms_ca = chain[p4].num_atoms_ca;
    memcpy(chain_local[3].pos,chain[p4].pos,sizeof(chain[p4].pos));

    debug("Computing local centroid.\n");
    chain_local[0].ComputeCentroid();
    chain_local[1].ComputeCentroid();
    chain_local[2].ComputeCentroid();
    chain_local[3].ComputeCentroid();

#ifndef NDEBUG
    debug("original:\n");
    chain[0].print_centroid();
    chain[1].print_centroid();
    chain[2].print_centroid();
    chain[3].print_centroid();
    // Compute the Centroids.
    debug("local:\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
    // exit(0);
#endif

    // for(int i=0; i<chain[p1].num_atoms_ca; i++){
    //     printf("%f ",chain_local[0].pos[i].x);
    // }
    // printf("\n");


    // Rotation. (locally)
    Vector t0,t1;
    Matrix R1,R2,R3;
    Matrix R21,R22,R23;
    // R4;
    // Matrix yx01,x12;
    double spot;


    Curvature curvleft;
    Curvature curvright;
    Curvature curvavg;


    // translate to origin
    t0 = translate_to_origin(chain_local,0,4);
    // printf("after local curvature translation:\n");
    // printf("local:\n");
    // chain_local[0].print_centroid();
    // chain_local[1].print_centroid();
    // chain_local[2].print_centroid();
    // chain_local[3].print_centroid();
    // printf("original:\n");
    // chain[0].print_centroid();
    // chain[1].print_centroid();
    // chain[2].print_centroid();
    // chain[3].print_centroid();
    // exit(0);


    // around Z to X axis, more X.
    // axisfrom Z(2) (no Z comp.), axisto X(0) (UNIT.X); sent to X.
    spot = chain_local[1].centroid.y;
    R1 = rotate_system_around_axis(chain_local,4,2,0,           \
                                   chain_local[0].centroid,     \
                                   chain_local[1].centroid,1);


#ifndef NDEBUG
    debug("pre-round 1.\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif

    // printf("y-spot vs. new y:   %f    %f\n",spot,chain_local[1].centroid.y);
    // Check1.
    if (fabs(chain_local[1].centroid.y) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 1 !!!\n");

        // TRANSPOSE
        R1.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R1,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R1,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();
        chain_local[3].ComputeCentroid();

        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        chain_local[3].print_centroid();
        debug("WARNING: had to exit at compute_curvature/1st rotation.\n");
#endif
        // exit(1);
    } // end of Check1.
    // exit(0);
    debug("----------Rotate #1 Done.-----------\n");



    spot = chain_local[1].centroid.z;
    R2 = rotate_system_around_axis(chain_local,4,1,0,           \
                                   chain_local[0].centroid,     \
                                   chain_local[1].centroid,-1);


#ifndef NDEBUG
    debug("pre-round 2.\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif

    // printf("y-spot vs. new y:   %f    %f\n",spot,chain_local[1].centroid.y);
    // Check1.
    if (fabs(chain_local[1].centroid.z) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 2 !!!\n");

        // TRANSPOSE
        R2.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R2,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R2,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();
        chain_local[3].ComputeCentroid();

        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        chain_local[3].print_centroid();
        debug("WARNING: had to exit at compute_curvature/2nd rotation.\n");
#endif
    } // end of Check1.
    debug("----------Rotate #2 Done.-----------\n");
    // exit(0);


    spot = chain_local[2].centroid.z;
    R3 = rotate_system_around_axis(chain_local,4,0,1,          \
                                    chain_local[1].centroid,    \
                                    chain_local[2].centroid,1);


#ifndef NDEBUG
    debug("pre-round 3.\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif


    // printf("z-spot vs. new z:   %f    %f\n",spot,chain_local[2].centroid.z);
    // Check1.
    if (fabs(chain_local[2].centroid.z) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 3 !!!\n");

        // TRANSPOSE
        R3.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R3,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R3,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();
        chain_local[3].ComputeCentroid();

        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        chain_local[3].print_centroid();
        debug("WARNING: had to exit at compute_curvature/3rd rotation.\n");
#endif
    } // end of Check1.
    debug("----------Rotate #3 Done.-----------\n");


#ifndef NDEBUG
    // print centroid.
    debug("curvature_check:\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif
    // quality check.
    printf("after_first_rotation(pos-2):\n");
    chain_local[2].print_centroid();

    curvleft.p1 = chain_local[0].centroid;
    curvleft.p2 = chain_local[1].centroid;
    curvleft.p3 = chain_local[2].centroid;
    curvleft.get_curvature();


    // // left (0-1-2) | p1,p2,p3?
    // if ((chain_local[0].centroid.z < 1.0) &&
    //     (chain_local[1].centroid.z < 1.0) &&
    //     (chain_local[2].centroid.z < 1.0)) {
    //     curvleft.p1 = chain_local[0].centroid;
    //     curvleft.p2 = chain_local[1].centroid;
    //     curvleft.p3 = chain_local[2].centroid;
    //     curvleft.get_curvature();
    // } else {
    //     printf("WARNING: bigz, out of place.\n");
    //     chain_local[0].print_centroid();
    //     chain_local[1].print_centroid();
    //     chain_local[2].print_centroid();
    //     curvleft.radius = 99999.0;
    //     curvleft.curvature = 99999.0;
    //     // exit(0);
    // }


    /* ---------------------------------------------------------
       Round 2.
       --------------------------------------------------------- */
    // translate to origin
    t1 = translate_to_origin(chain_local,1,4);
    // printf("after local curvature translation:\n");

#ifndef NDEBUG
    debug("pre-round 2-1.\n");
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
    debug("original:\n");
    chain[0].print_centroid();
    chain[1].print_centroid();
    chain[2].print_centroid();
    chain[3].print_centroid();
    // exit(0);
#endif

    // around Z to X axis, more X.
    // axisfrom Z(2) (no Z comp.), axisto X(0) (UNIT.X); sent to X.
    spot = chain_local[2].centroid.y;
    R21 = rotate_system_around_axis(chain_local,4,2,0,           \
                                    chain_local[1].centroid,     \
                                    chain_local[2].centroid,1);

#ifndef NDEBUG
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif

    // printf("y-spot vs. new y:   %f    %f\n",spot,chain_local[1].centroid.y);
    // Check1.
    if (fabs(chain_local[2].centroid.y) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 2-1 !!!\n");

        // TRANSPOSE
        R21.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R21,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }
        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R21,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();
        chain_local[3].ComputeCentroid();

        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        chain_local[3].print_centroid();
        debug("WARNING: had to exit at compute_curvature/21 rotation.\n");
#endif
        // exit(1);
    } // end of Check1.
    debug("----------Rotate #2-1 Done.-----------\n");
    // exit(0);

    spot = chain_local[2].centroid.z;
    R22 = rotate_system_around_axis(chain_local,4,1,0,          \
                                    chain_local[1].centroid,        \
                                    chain_local[2].centroid,-1);

#ifndef NDEBUG
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif

    // printf("y-spot vs. new y:   %f    %f\n",spot,chain_local[1].centroid.y);
    // Check1.
    if (fabs(chain_local[2].centroid.z) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 2-2 !!!\n");

        // TRANSPOSE
        R22.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R22,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }
        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R22,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();
        chain_local[3].ComputeCentroid();
        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        chain_local[3].print_centroid();
        debug("WARNING: had to exit at compute_curvature/2nd rotation.\n");
#endif
    } // end of Check1.
    debug("----------Rotate # 2-2 Done.-----------\n");


    spot = chain_local[3].centroid.z;
    R23 = rotate_system_around_axis(chain_local,4,0,1,           \
                                    chain_local[2].centroid,     \
                                    chain_local[3].centroid,1);

#ifndef NDEBUG
    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif

    // printf("z-spot vs. new z:   %f    %f\n",spot,chain_local[2].centroid.z);
    // Check1.
    if (fabs(chain_local[3].centroid.z) > fabs(spot)) {
        // need to double rotate back.
        debug("Check 2-3 !!!\n");

        // TRANSPOSE
        R23.transpose();

        Vector newpos;
        // rotate once.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R23,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }
        // rotate twice.
        for(int i=0; i<4; i++){
            for(int j=0; j<chain_local[i].num_atoms_ca; j++){
                // printf("%d",j);
                newpos = prod_matrix_vector(R23,chain_local[i].pos[j]);
                chain_local[i].pos[j].x = newpos.x;
                chain_local[i].pos[j].y = newpos.y;
                chain_local[i].pos[j].z = newpos.z;
            }
            // printf("\n");
        }

#ifndef NDEBUG
        // compute centroid.
        chain_local[0].ComputeCentroid();
        chain_local[1].ComputeCentroid();
        chain_local[2].ComputeCentroid();
        chain_local[3].ComputeCentroid();

        // print centroid.
        chain_local[0].print_centroid();
        chain_local[1].print_centroid();
        chain_local[2].print_centroid();
        chain_local[3].print_centroid();
        debug("WARNING: had to exit at compute_curvature/3rd rotation.\n");
#endif
    } // end of Check1.
    debug("----------Rotate # 2-3 Done.-----------\n");

#ifndef NDEBUG
    // print centroid.
    debug("curvature_check:\n");

    chain_local[0].print_centroid();
    chain_local[1].print_centroid();
    chain_local[2].print_centroid();
    chain_local[3].print_centroid();
#endif
    // quality check.
    printf("after_second_rotation(pos-3):\n");
    chain_local[3].print_centroid();

    curvright.p1 = chain_local[1].centroid;
    curvright.p2 = chain_local[2].centroid;
    curvright.p3 = chain_local[3].centroid;
    curvright.get_curvature();

    // // right (1-2-3)
    // if ((chain_local[1].centroid.z < 11.0) &&
    //     (chain_local[2].centroid.z < 11.0) &&
    //     (chain_local[3].centroid.z < 11.0)) {
    // } else {
    //     printf("WARNING: bigz, out of place.\n");
    //     chain_local[1].print_centroid();
    //     chain_local[2].print_centroid();
    //     chain_local[3].print_centroid();
    //     curvright.radius = 99999.0;
    //     curvright.curvature = 99999.0;
    //     // exit(1);
    // }


    // spot = chain_local[3].centroid.z;
    // R4 = rotate_system_around_axis(chain_local,4,0,1,          \
    //                                 chain_local[2].centroid,    \
    //                                 chain_local[3].centroid,1);

    // chain_local[0].print_centroid();
    // chain_local[1].print_centroid();
    // chain_local[2].print_centroid();
    // chain_local[3].print_centroid();

    // // printf("z-spot vs. new z:   %f    %f\n",spot,chain_local[3].centroid.z);
    // // Check1.
    // if (fabs(chain_local[3].centroid.z) > fabs(spot)) {
    //     // need to double rotate back.
    //     printf("Check4 !!!\n");

    //     // TRANSPOSE
    //     R4.transpose();

    //     Vector newpos;
    //     // rotate once.
    //     for(int i=0; i<4; i++){
    //         for(int j=0; j<chain_local[i].num_atoms_ca; j++){
    //             // printf("%d",j);
    //             newpos = prod_matrix_vector(R4,chain_local[i].pos[j]);
    //             chain_local[i].pos[j].x = newpos.x;
    //             chain_local[i].pos[j].y = newpos.y;
    //             chain_local[i].pos[j].z = newpos.z;
    //         }
    //         // printf("\n");
    //     }

    //     // rotate twice.
    //     for(int i=0; i<4; i++){
    //         for(int j=0; j<chain_local[i].num_atoms_ca; j++){
    //             // printf("%d",j);
    //             newpos = prod_matrix_vector(R4,chain_local[i].pos[j]);
    //             chain_local[i].pos[j].x = newpos.x;
    //             chain_local[i].pos[j].y = newpos.y;
    //             chain_local[i].pos[j].z = newpos.z;
    //         }
    //         // printf("\n");
    //     }

    //     // compute centroid.
    //     chain_local[0].ComputeCentroid();
    //     chain_local[1].ComputeCentroid();
    //     chain_local[2].ComputeCentroid();
    //     chain_local[3].ComputeCentroid();

    //     // print centroid.
    //     chain_local[0].print_centroid();
    //     chain_local[1].print_centroid();
    //     chain_local[2].print_centroid();
    //     chain_local[3].print_centroid();

    //     printf("WARNING: had to exit at compute_curvature/4th rotation.\n");
    //     // exit(1);
    // } // end of Check1.
    // printf("----------Rotate #4 Done.-----------\n");
    // // exit(0);



    // Vector onefour;
    // onefour = get_vector(chain_local[0].centroid,chain_local[3].centroid);
    // printf("\nthe_onefour_vector:");
    // onefour.print_Vector();

    // Rz.print_Matrix();
    // yx01.print_Matrix();
    // x12.print_Matrix();

    // printf("after Z to X:\n");
    // printf("after Y to X:\n");
    // printf("local:\n");
    // chain_local[0].print_centroid();
    // chain_local[1].print_centroid();
    // chain_local[2].print_centroid();
    // chain_local[3].print_centroid();
    // printf("original:\n");
    // chain[0].print_centroid();
    // chain[1].print_centroid();
    // chain[2].print_centroid();
    // chain[3].print_centroid();
    // exit(0);


    // around X to Y axis, less Z.
    // axisfrom ()
    // Rx = rotate_system_around_axis(chain_local,4,0,1,       \
    //                                chain_local[0].centroid,     \
    //                                chain_local[3].centroid,1);
    // Rz.print_Matrix();
    // printf("after X to Y:\n");
    // printf("local:\n");
    // chain_local[0].print_centroid();
    // chain_local[1].print_centroid();
    // chain_local[2].print_centroid();
    // chain_local[3].print_centroid();
    // printf("original:\n");
    // chain[0].print_centroid();
    // chain[1].print_centroid();
    // chain[2].print_centroid();
    // chain[3].print_centroid();
    // exit(0);


    // // print Centroids
    // Vector proto_midpoint,proto_centroid;
    // proto_centroid = get_centroid_from_centroids(chain_local,0,3);
    // proto_midpoint = midpoint2(chain_local[0].centroid,chain_local[4-1].centroid);
    // proto_centroid.print_Vector();
    // proto_midpoint.print_Vector();
    // // around Y to
    // Ry = rotate_system_around_axis(chain_local,4,1,2,            \
    //                                proto_midpoint,               \
    //                                proto_centroid,1);
    // // exit(0);

    // printf("after X to Y:\n");
    // chain[0].print_centroid();
    // chain[1].print_centroid();
    // chain[2].print_centroid();
    // chain[3].print_centroid();
    // exit(0);

    // print before
#ifndef NDEBUG
    printf("curv left and right, before:\n");
    curvleft.print_rc();
    curvright.print_rc();
#endif


    // get average
    if ((curvleft.radius < 2000.0) && (curvright.radius < 2000.0)) {
        curvavg.radius = (curvleft.radius + curvright.radius) * 0.5;
    } else if ((curvleft.radius > 2000.0) || (curvright.radius > 2000.0)) {
        curvavg.radius = 2000.01;
    }

    if ((curvleft.curvature < 500.0) && (curvright.curvature < 500.0)) {
        curvavg.curvature = (curvleft.curvature + curvright.curvature) * 0.5;
    } else if ((curvleft.curvature > 500.0) || curvright.curvature > 500.0) {
        curvavg.curvature = 500.01;
    }



#ifndef NDEBUG
    debug("average:\n");
    curvavg.print_rc();
#endif

#ifndef NDEBUG
    debug("difference:\n");
    debug("\tradius: %f \n\tcurvature: %f\n",curvleft.radius - curvavg.radius,curvleft.curvature - curvavg.curvature);
#endif


    // int change;
    // change = 0;

    // // radius > 999.0 -> 999.0
    // if (curvavg.radius > 999.0){
    //     curvavg.radius = 999.0;
    //     change += 1;
    // }

    // // curvature: nan -> 0
    // if (isnan(curvavg.curvature) == 1){
    //     curvavg.curvature = 0.0;
    //     change += 1;
    // }

    // if (change > 0) {
    //     curvavg.print_rc();
    //     printf("modify: averages' limits were hit!!!!!  WARNING!!\n");
    // } else {
    //     printf("modify: actual average used!\n");
    // }



    // FPRINTF
    // fprintf(fp_curvature_local,"%8.2f ",local_curv.curvature);
    // fprintf(fp_radius_of_curvature_local,"%8.2f ",local_curv.radius);

    // Check that p1,p2,p3 are in curvleft,curvright
    // curvleft.print_points();
    // curvright.print_points();
    // exit(0);

    fprintf(fp_curv,"%6.2f ",curvavg.curvature);
    fprintf(fp_curv,"%6.3f ",curvleft.p3.z);
    fprintf(fp_curv,"%6.3f   ",curvright.p3.z);

    fprintf(fp_rad,"%6.2f ",curvavg.radius);
    fprintf(fp_rad,"%6.3f ",curvleft.p3.z);
    fprintf(fp_rad,"%6.3f   ",curvright.p3.z);


    // delete
    delete [] chain_local;
    return curvavg;

    // Vector midhighlow;
    // for(int i=1; i<chains_to_use-1; i+=2) {
    //     printf("%d\n",i);
    //     midhighlow = midpoint2(chain_later[i].pos[chain_later[i].pf_pos_high],\
    //                            chain_later[i+1].pos[chain_later[i].pf_pos_low]);
    //     // printf("%f %f %f\n",midhighlow.x,midhighlow.y,midhighlow.z);
    //     radius_curv,degrees = curvature(chain_later[i].centroid,
    //                              midhighlow,
    //                              chain_later[i+1].centroid);
    //     fprintf(fp_curvature_local,"%8.2f ",degrees);
    //     fprintf(fp_radius_of_curvature_local,"%8.2f ",radius_curv);
    // }

    // local_curv.print_points();
    // for(int i=1; i<chains_to_use-1; i+=2) {
        // printf("\nlocal_curvature: %d\n",i);
        // curvature(chain,i) chain=entire system, i=beta monomer
        // curvatureCH(chain_later,i,&radius_curv,&degrees);
        // fprintf(fp_curvature_local,"%8.2f ",degrees);
        // fprintf(fp_radius_of_curvature_local,"%8.2f ",radius_curv);
    // }

}

Vector midpoint_of_indices_of_2_chains(Chain *chain,int c1,int c2,int &idx1, int &idx2){
    // c1: chain 1.
    // c2: chain 2.
    debug("inside closest indices!\n");
    debug("chains: %d %d\n",c1,c2);
    double dist,dist_lowest;
    dist_lowest = distance(chain[c1].pos[0],chain[c2].pos[0]);

    for(int i=0; i<chain[c1].num_atoms_ca; i++){
        for(int j=0; j<chain[c2].num_atoms_ca; j++){
            dist = distance(chain[c1].pos[i],chain[c2].pos[j]);
            if(dist < dist_lowest){
                dist_lowest = dist;
                idx1 = i;
                idx2 = j;
            }
        }
    }

    // std::vector<int> list5;
    // double dist,dist_in1;
    // int count = 0;
    // // std::vector<int> list3;
    // // list3.front() = 0;

    // for(int i=0; i<chain[c1].num_atoms_ca; i++){
    //     printf("%d\n",i);
    //     dist = distance(chain[c1].pos[i],chain[c2].centroid);
    //     // if(list5.size()<=4){
    //     //     list5.push_back(i);
    //     //     continue;
    //     // }

    //     if(list5.size()>4){
    //         printf("at 5..\n");
    //         exit(0);
    //         for(int j=0; j<list5.size(); j++){
    //             dist_in1 = distance(chain[c1].pos[list5[j]],chain[c2].centroid);
    //             if(dist < dist_in1){
    //                 printf("current_dist: %f\n",dist_in1);
    //                 printf("lower_new_dist: %f\n",dist);
    //                 list5[count] = i;
    //                 count += 1;

    //                 break;
    //             }
    //         }
    //         // break;
    //     }
    // }


    // debug("the closest were %d to %d\n",idx1,idx2);
    debug("chains: %d %d, >>the_index: %d  to  %d, at dist: %f\n", \
          c1,c2,chain[c1].indices[idx1],chain[c2].indices[idx2],\
          dist_lowest);

    Vector mpoint;
    mpoint = midpoint2(chain[c1].pos[idx1],chain[c2].pos[idx2]);

    debug("%f %f %f\n",mpoint.x,mpoint.y,mpoint.z);
    return mpoint;
}
int get_pos_closest_to_point(Chain *chain,int c,Vector P){

    int index;
    double dist,dist_low;

    dist_low = distance(chain[c].pos[0],P);
    for(int i=0; i<chain[c].num_atoms_ca; i++){
        dist = distance(chain[c].pos[i],P);
        if(dist<dist_low){
            dist_low = dist;
            index = i;
        }
    }

    debug(" to point P: %f %f %f\n",P.x,P.y,P.z);
    debug("closest POS is: %d\n",index);
    debug("closest INDEX is: %d\n",chain[c].indices[index]);
    return index;
}
Vector translate_to_origin(Chain *chain,int num_chain,int max_num){
    //
    Vector newpos;
    newpos.x = chain[num_chain].centroid.x + 0.001;
    newpos.y = chain[num_chain].centroid.y + 0.001;
    newpos.z = chain[num_chain].centroid.z + 0.001;

    // rotate
    for(int i=0; i<max_num; i++){
        for(int j=0; j<chain[i].num_atoms_ca; j++){
            // printf("%d",j);
            // newpos = prod_matrix_vector(rotZ,chain[i].pos[j]);
            chain[i].pos[j].x -= newpos.x;
            chain[i].pos[j].y -= newpos.y;
            chain[i].pos[j].z -= newpos.z;
        }
        // printf("\n");
    }


    for ( int i=0; i<max_num; i++) {
        chain[i].ComputeCentroid();
    }

    return newpos;
}
Matrix rotate_system_around_axis(Chain *chain,int max_num_chains,int axisfrom,int axisto,Vector O,Vector E,int d){
    // axis:  x,0;  y,1;  z,2;
    // from z to y,  2 to 1

    Vector V,norm_V,UNIT;
    get_vector(O,E,&V);

    if(axisfrom == 0){
        V.x = 0.0;
    } else if (axisfrom == 1){
        V.y = 0.0;
    } else if (axisfrom == 2){
        V.z = 0.0;
    }
    normalize(V,&norm_V);

    // unit vector
    if(axisto == 1){
        // Y
        UNIT.x = 0.0;
        UNIT.y = 1.0;
        UNIT.z = 0.0;
    } else if (axisto == 2){
        // Z
        UNIT.x = 0.0;
        UNIT.y = 0.0;
        UNIT.z = 1.0;
    } else {
        // X (default)
        UNIT.x = 1.0;
        UNIT.y = 0.0;
        UNIT.z = 0.0;
    }

    double costheta,theta; // ,sintheta; // with Y
    costheta = get_costheta(norm_V,UNIT);

    // sintheta = get_sintheta(norm_V,Y);
    debug("costheta: %2.3f\n",costheta);
    // printf("[costheta,sintheta]: %2.3f  %2.3f\n",costheta,sintheta);


    theta = acos(costheta);
    debug("theta: %f\n",theta);
    Matrix Rot;

    // direction
    double D;
    D = (double)d;

    debug("double: __%d__\n",d);
    debug("double: __%f__\n",D);

    // exit(0);
    Rot.build_rotation(theta * -1.0 * D,axisfrom); // X
#ifndef NDEBUG
    Rot.print_Matrix();
#endif


    Vector newpos;
    // rotate
    for(int i=0; i<max_num_chains; i++){
        for(int j=0; j<chain[i].num_atoms_ca; j++){
            // printf("%d",j);
            newpos = prod_matrix_vector(Rot,chain[i].pos[j]);
            chain[i].pos[j].x = newpos.x;
            chain[i].pos[j].y = newpos.y;
            chain[i].pos[j].z = newpos.z;
        }
        // printf("\n");
    }



    for ( int i=0; i<max_num_chains; i++) {
        chain[i].ComputeCentroid();
    }

    return Rot;
}

void get_xy_plane(Chain *chain,int max_num_chains,int num_chain){
    // http://www.maplesoft.com/support/help/AddOns/view.aspx?path=MathApps/EquationofaPlane3Points
    // chain indices: like 0,1,2 chains in system, or 2,3,4
    int a,b,c;
    a = num_chain;
    b = num_chain + 1;
    c = num_chain + 2;
    printf("a:%d  b:%d  c:%d\n",a,b,c);
    printf("%f   %f   %f\n",chain[a].centroid.x,chain[a].centroid.y,chain[a].centroid.z);
    printf("%f   %f   %f\n",chain[b].centroid.x,chain[b].centroid.y,chain[b].centroid.z);
    printf("%f   %f   %f\n",chain[c].centroid.x,chain[c].centroid.y,chain[c].centroid.z);

    // ab is b-a vector
    // ac is c-a vector
    // V is the cross product of ab,ac. ( ab X ac )
    Vector ab,ac,V;
    Vector norm_V;

    // Vector prime_a, prime_b, prime_c;
    // http://stackoverflow.com/questions/6264664/transform-3d-points-to-d2
    // A' = A / ||V||
    // B' = B / ||V||
    // C' = C / ||V||
    // ||V|| = (A2+B2+C2)1/2
    double prime_a, prime_b, prime_c;


    // get ab,ac
    // void get_vector ( Vector v1, Vector v2, Vector *v3 ) {
    get_vector(chain[a].centroid,chain[b].centroid,&ab); // ab
    get_vector(chain[a].centroid,chain[c].centroid,&ac); // ac
    V = cross_product(ab,ac);

    // Equation of a plane. Ax + By +Cz + d = 0.
    // float d;
    // d = (V.x * chain[a].centroid.x \
    //     + V.y * chain[a].centroid.y \
    //      + V.z * chain[a].centroid.z) * -1.0;
    // printf("A:%f   B:%f   C:%f   d:%f\n",V.x,V.y,V.z,d);

    // normalize & magnitude V
    double magnitude_V;
    double theta,costheta,sintheta,m1costheta;
    magnitude_V = magnitude(V);
    prime_a = V.x / magnitude_V;
    prime_b = V.y / magnitude_V;
    prime_c = V.z / magnitude_V;
    printf("primes: %f   %f   %f\n",prime_a,prime_b,prime_c);

    printf("c_prime: %f\n",prime_c);
    theta = acos(prime_c);
    sintheta = sin(theta);
    costheta = prime_c;
    printf("theta: %f\n",theta);

    Vector Vp, Z, R, norm_R;
    Vp.x = prime_a;
    Vp.y = prime_b;
    Vp.z = prime_c;

    Z.x = 0.0;
    Z.y = 0.0;
    Z.z = 1.0;

    R = cross_product(Vp,Z);
    printf("R:  %f   %f   %f\n",R.x,R.y,R.z);

    normalize(R,&norm_R);
    printf("norm_R:  %f   %f   %f\n",norm_R.x,norm_R.y,norm_R.z);


    // https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
    // Rotation Vector:
    // Theta = theta . e   // (e=Vp)

    // Rotating a vector. (Rodrigues' rotation formula)
    // v will be the points ???
    // V_rot = (cos theta)v + (sintheta)(e X v) + (1 - costheta)(e . v) e
    // vectors:  v_rot, v_ctv, v_sev,v_eve;
    Vector v_rot, v_ctv, v_sev, v_eve;


}

// testing
// std::vector<int> get_protofilaments(Chain *chain,int chains_to_use) {
std::vector<int> get_monomers_with_no_southern_neighbor(Chain *chain,int chains_to_use) {
    printf("getting protofilaments .. from the %d chains.\n",chains_to_use);


    // std::vector<int> pf (1,-1); // start with a list of length 1, initialized to -1.
    std::vector<int> pf; // empty vector.
    std::vector<int>::iterator it;

    for( int i = 0; i<chains_to_use; i++ ) {

        // printf("NESW_(cw): %d %d %d %d\n",chain[i].LonN,chain[i].LatE,\
        //        chain[i].LonS,chain[i].LatW);

        if (chain[i].LonS == -1) {

            // printf("%d\n",i);
            // printf("\tNESW_(cw): %d %d %d %d\n",chain[i].LonN,chain[i].LatE,\
            //        chain[i].LonS,chain[i].LatW);

            // pf.insert(pf.begin(),i); // becomes 24 22 20 18 ...
            pf.insert(pf.end(),i);
            // pf.insert

        } // if
    } // for (chains_to_use)

    // std::cout << "protofilaments identified:" << std::endl; // has an implicit \n
    printf("monomers with no southern neighbor identified:\n");
    for (it=pf.begin(); it<pf.end(); it++) {
        // std::cout << *it << ' ' << std::endl;
        printf("%d ",*it);
    }
    std::cout << ' ' << std::endl; // has the implicit \n

    return pf;
}

std::vector<int> get_protofilament(Chain *chain,int pf) {
    printf("getting protofilament .. from the %d chain.\n",pf);

    int northern_neighbor=-2;
    int current_chain = pf;

    // std::vector<int> pf (1,-1); // start with a list of length 1, initialized to -1.
    std::vector<int> proto; // empty vector.
    std::vector<int>::iterator it;

    // for( int i = 0; i<chains_to_use; i++ ) {

    // }

    proto.insert(proto.end(),pf);

    while (northern_neighbor != -1) {
        northern_neighbor = chain[current_chain].LonN;

        if (northern_neighbor != -1) {
            proto.insert(proto.end(),northern_neighbor);
        }
        printf("chain %d is south of %d.\n",current_chain,northern_neighbor);
        current_chain = northern_neighbor;
    }

    return proto;
}
