// chi.h
#ifndef _CHI_H_
#define _CHI_H_

// libraries:
/* #include <stdio.h> */

// Definitions:
/* #define BUFFERSIZE 900 */

/* ---------------------------------------------------------
   Class declarations
   --------------------------------------------------------- */
/* #ifndef CHAIN_H */
/* #define CHAIN_H */
#include "chain.h"
/* #endif */


/* ---------------------------------------------------------
   function declarations
   --------------------------------------------------------- */
double get_chi_global(Chain *chain_ref,Chain *chain_later,int L1,int L2);
void get_chi_by_residue(Chain *chain_ref,Chain *chain_later,int L1,int L2,FILE *fp);

/* void get_xy_plane(Chain *chain,int max_num_chains,int num_chain); */
/* double curvature(Vector p1,Vector p2,Vector p3); */
/* void translate_to_origin(Chain *chain,int num_chain,int max_num); */
/* void rotate_system_around_z_to_y(Chain *chain,int begin,int end); */
/* void rotate_system_around_axis(Chain *chain,int max_num_chains,\ */
/*                                int axisfrom,int axisto,Vector Origin,Vector Endpoint); */
/* Vector midpoint_of_indices_of_2_chains(Chain *chain,int c1,int c2,int &idx1, int &idx2); */
/* int get_pos_closest_to_point(Chain *chain,int c,Vector v); */
#endif
