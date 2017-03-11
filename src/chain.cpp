// chain.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
// #include <string>
#include <string> // strcpy, memcpy
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

/* ---------------------------------------------------------
   other headers
   --------------------------------------------------------- */
#include <chain.h>
#include <debug.h>
#include "boost/format.hpp"

// Chain chain;
// #include "include/molecule.h"
Molecule mol;

#include "md.h"
Vector pos;


// #include "dcdio.h"

/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
void build_contact_map( Molecule *mol ) {
    debug("building contact map for %s\n",mol->filename);

    int i;
    int j; /* -> 40, MX_CONTACTS */

    /* printf */
    for ( i = 0; i < mol->num_atoms_ca; i++ ) {
        /* printf("%d\n",i); */
        for ( j = 0; j < MAX_CONTACTS; j ++ ) {
            /* if ( mol->contacts[i][j].cresid != 0 || mol->contacts[i][j].distance != 0.0) { */
            /* printf("%d %f\n",mol->contacts[i][j].cresid,mol->contacts[i][j].distance); */
            mol->contacts[i][j].cresid = -1;
            mol->contacts[i][j].distance = 0.0;
            /* } */
        }
    }

    /* for ( j = 0; j < 40; j ++ ) { */
    /* } */

    /* for ( i = 0; i < mol->num_atoms_ca; i++ ) { */
    /* if ( i) */
    /* printf("%d\n",i); */
    /* for ( j = 0; j < 40; j ++ ) { */
    /* } */
    /* } */

    /* int i; // for loop iterator */
    /* int j; // for loop iterator (2) */
    int count; // while loop
    double dist;

    for( i = 0; i < mol->num_atoms_ca; i++ ) {
        for( j = i; j < mol->num_atoms_ca; j++) {

            if ( i >= j-2 && i <= j+2) {
                // exclude j = i +/- 2
                /* printf("%d\t no --> %d j: %d %d\n",i,j,j-2,j+2); */
            }
            else {
                /* printf("%d\n",j); */
                dist = distance( mol->pos[i], mol->pos[j]);
                /* printf("%8.4f\n",dist); */

                if ( dist <= 8.0 ) {
                    /* printf("criterion-> %d %d %8.4f\n",i,j,dist); */


                    /* SINGLE */
                    /* for( count = 0; count < 40; count++ ) { */
                    /*     if ( mol->contact_map[i][count].cresid != -1 ) { */
                    /*         continue; */
                    /*     } else { */
                    /*         mol->contact_map[i][count].cresid = j; */
                    /*         mol->contact_map[i][count].distance = dist; */
                    /*         break; */
                    /*     } */
                    /* } */


                    /* DOUBLE */
                    for( count = 0; count < MAX_CONTACTS; count++ ) {

                        if (mol->contacts[i][count].cresid != -1 ) {
                            continue;
                        } else {
                            mol->contacts[i][count].cresid = j;
                            mol->contacts[i][count].distance = dist;
                            // mol->contacts[i][k].cresid = j;
                            // mol->contacts[i][k].distance = dist;
                            mol->contacts[i][count].index = i + mol->index;
                            mol->contacts[i][count].cindex = j + mol->index;

                            break;
                        }
                    }
                    /* and j */
                    for( count = 0; count < MAX_CONTACTS; count++ ) {

                        if (mol->contacts[j][count].cresid != -1 ) {
                            continue;
                        } else {
                            mol->contacts[j][count].cresid = i;
                            mol->contacts[j][count].distance = dist;
                            // mol->contacts[j][l].cresid = i;
                            // mol->contacts[j][l].distance = dist;
                            mol->contacts[j][count].index = j + mol->index;
                            mol->contacts[j][count].cindex = i + mol->index;

                            break;
                        }
                    }
                    /* break; */
                    /* mol->topology[i][ */
                    /* } else { */
                    /*     count */
                    /* } */
                    /* } */
                } /* END 8.0 */
            }
        }
    }
    /* printf("the count was: %d\n", count); */

    int tally;
    int tally_total = 0;

    for ( i = 0; i < mol->num_atoms_ca; i++ ) {
        /* printf("%d\n",i); */
        tally = 0;
        for ( j = 0; j < MAX_CONTACTS; j ++ ) {
            /* if ( mol->contacts[i][j].cresid != 0 || mol->contacts[i][j].distance != 0.0) { */
            if ( mol->contacts[i][j].cresid != -1) {
                /* printf("%d %f %d\n",mol->contacts[i][j].cresid,mol->contacts[i][j].distance,j); */
                tally++;
            }
        }

        /* PRINTHERE */
        /* printf("tally for %d: %d\n",i,tally); */
        mol->tally_num_contacts_2x[i] = tally;
        tally_total += tally;
    }
    debug("tally_total: %d \t %d\n",tally_total,tally_total/2);
    mol->total_num_contacts_2x = tally_total;
}

void fprintf_count_of_contacts_2x_persist_from_reference( Molecule *m1, Molecule *m2 ) {
    // Molecule m1: reference
    // Molecule m2: contacts lost
    int i; //
    int j; //
    int k; //
    int count = 0; // contacts, for comparison
    /* int len_array_num_contacts = sizeof(m1->num_contacts) / sizeof(m1->num_contacts[0]); */
    /* int len_array_num_contacts_2x = sizeof(m1->num_contacts_2x) / sizeof(m1->num_contacts_2x[0]); */

    printf("generating: ref_contact_map.dat....\n");
    // uncomment 7/8 printf
    // contact_map ==? contact_map_2x
    // num_contacts == ? num_contacts_2x



    for ( i = 0; i < m1->num_atoms_ca; i++ ) {
        count = 0;

        /* printf("%d:",i); */
        /* printf("\nref:"); */


        /* PRINT HERE */
        /* for ( j = 0; j < m1->tally_num_contacts_2x[i]; j ++ ) { */
        /*     /\* printf("%d",j); *\/ */
        /*     printf(" %d-%d ",i,m1->contacts[i][j].cresid); */
        /*     /\* printf("%d"); *\/ */
        /* } */


        /* printf("\n"); */
        /* printf("com:"); */


        for ( j = 0; j < m2->tally_num_contacts_2x[i]; j ++ ) {

            /* printf(" %d-%d",i,m2->contact_map_2x[i][j].cresid); */

            for ( k = 0; k < m1->tally_num_contacts_2x[i]; k++ ) {

                /* PRINT HERE */
                /* printf("%d\n",k); */
                if ( m2->contacts[i][j].cresid == m1->contacts[i][k].cresid ) {
                    count++;

                    /* PRINT HERE */
                    /* printf("%d %d %d\n",i,m1->contacts[i][k].cresid,m2->contacts[i][j].cresid); */
                }

            } /* end for */
        } /* end for */

        /* printf("\n"); */
        printf("%d",count);
        m2->num_contacts_2x_persist_ref[i] = count;
    }

    int total_persisting_contacts = 0;
    for ( i = 0; i < m2->num_atoms_ca; i++ ) {
        total_persisting_contacts += m2->num_contacts_2x_persist_ref[i];
    }
    printf("\nPERSISTING_contact_total: %s --> %s was: %d, %d\n", \
           m1->filename,m2->filename,total_persisting_contacts,total_persisting_contacts/2);

    // PRINTF
    /* for ( i = 0; i < m2->num_atoms_ca; i++ ) { */
    /*     printf("%4d",m2->num_contacts_2x_persist_ref[i]); */
    /* } */

    FILE * fp_new_contacts;
    fp_new_contacts = fopen ("ref_contact_map_C.dat", "a+");

    // FPRINTF
    /* for ( i = 0; i < m2->num_atoms_ca; i++ ) { */
    for ( i = 0; i < m2->num_atoms_ca; i++ ) {
        /* printf(" %4d",m2->num_contacts_2x_persist_ref[i]); */

        /* Main */
        fprintf(fp_new_contacts," %4d",m2->num_contacts_2x_persist_ref[i]);

        /* temp */
        /* fprintf(fp_new_contacts,"%2d\n",m2->num_contacts_2x_persist_ref[i]); */
    }
    /* } */
    fprintf(fp_new_contacts,"\n");
    fclose(fp_new_contacts);
}

void fprintf_bond_vector_angle_with_tension_vector( Molecule *m1, Molecule *m2 ) {
    // Molecule m1: reference
    // Molecule m2: with tension


    /* the FENE potential */
#define R_limit  2.0  /*8.0*/     /* this is the R0 parameter in FENE model */
#define kspring_cov /*40.0*/ 20.0 /* spring constant in FENE model (for the covalently bound beads) */
#define prefactor 70.0



    /* end to end in m2 (the changing molecule) */
    Vector endtoend;
    Vector endtoend_norm;
    double magnitude = 0.0;

    endtoend.x = m2->pos[m2->num_atoms_ca-1].x - m2->pos[0].x;
    endtoend.y = m2->pos[m2->num_atoms_ca-1].y - m2->pos[0].y;
    endtoend.z = m2->pos[m2->num_atoms_ca-1].z - m2->pos[0].z;

    printf("position_0: %f %f %f\n",m2->pos[0].x,m2->pos[0].y,m2->pos[0].z);
    printf("position_last: %f %f %f\n",m2->pos[m2->num_atoms_ca-1].x,m2->pos[m2->num_atoms_ca-1].y,m2->pos[m2->num_atoms_ca-1].z);
    // magnitude = normalize(endtoend,&endtoend_norm);
    normalize(endtoend,&endtoend_norm);
    // printf("endtoend_magnitude: %f\n",magnitude);
    printf("endtoend_vector: %f %f %f\n",endtoend.x,endtoend.y,endtoend.z);

    printf("computing: bondvector angle with tension vector\n");
    printf("computing: tension using R_limit(%f) and kspring_cov(%f)\n", R_limit,kspring_cov);


    int i;

    /* costheta */
    double costheta = 0.0;
    double arr_costheta[m2->num_atoms_ca];

    /* For 0--1--2 */
    /* Vector bond10; /\* 0 minus 1 *\/ */
    /* Vector bond12; /\* 2 minus 1 *\/ */
    Vector bond01; /* 1 minus 0 */
    Vector bond01_norm;

    /* for ( i = 0; i < m2->num_atoms_ca-1; i++ ) { */
    for ( i = 0; i < m2->num_atoms_ca; i++ ) {

        if ( i == 0 || i == m2->num_atoms_ca ) {
            continue;
        } else {
            bond01.x = m2->pos[i].x - m2->pos[i-1].x;
            bond01.y = m2->pos[i].y - m2->pos[i-1].y;
            bond01.z = m2->pos[i].z - m2->pos[i-1].z;
        }

        normalize(bond01,&bond01_norm);

        /* compute cos theta */
        costheta = get_costheta(bond01_norm,endtoend_norm);
        /* PRINT HERE */
        /* printf("costheta: %f\n",costheta); */
        arr_costheta[i] = costheta;
    }


    FILE * fp_bond_vector_angle_with_tension_vector;
    fp_bond_vector_angle_with_tension_vector = fopen ("costheta_bond_vector_with_tension_vector.dat", "a+");

    // FPRINTF
    fprintf(fp_bond_vector_angle_with_tension_vector,"# endtoend_vector: %f %f %f\n",endtoend.x,endtoend.y,endtoend.z);
    fprintf(fp_bond_vector_angle_with_tension_vector,"# norm_endtoend: %f %f %f\n",endtoend_norm.x,endtoend_norm.y,endtoend_norm.z);

    for ( i = 0; i < m2->num_atoms_ca-1; i++ ) {
        /* Main */
        fprintf(fp_bond_vector_angle_with_tension_vector,"%7.4f ",arr_costheta[i]);
    }
    fprintf(fp_bond_vector_angle_with_tension_vector,"\n");
    fclose(fp_bond_vector_angle_with_tension_vector);




    /* tension */
    double bond_length_ref[m1->num_atoms_ca - 1];
    double bond_length[m2->num_atoms_ca - 1];
    double diff_bond_length[m2->num_atoms_ca - 1];

    /* /\* UNUSED POTENTIAL: 1 of 2. *\/ */
    /* double potential[m2->num_atoms_ca - 1]; */
    double force[m2->num_atoms_ca - 1];
    double tforce[m2->num_atoms_ca -1]; /* force component in line with tension (endtoend); */


    /* tension */
    for ( i = 0; i < m1->num_atoms_ca - 1; i++ ) {
    /* for 9 atoms, 8 bonds, i -> 0 -- 7 */
        bond_length_ref[i] = distance(m1->pos[i],m1->pos[i+1]);
        bond_length[i] = distance(m2->pos[i],m2->pos[i+1]);
        diff_bond_length[i] = bond_length[i] - bond_length_ref[i];


        /* /\* UNUSED POTENTIAL: 2 of 2. *\/ */
        /* /\* pot = - (kspring_cov/2.0) * (R_limit*R_limit) * log(1. - dR*dR/(R_limit*R_limit)); *\/ */
        /* potential[i] = - (kspring_cov/2.0) * (R_limit*R_limit) *        \ */
        /*     log(1.0 - diff_bond_length[i] * diff_bond_length[i] / (R_limit * R_limit)); */

        /* force = kspring_cov*dR/(1. - dR*dR/(R_limit*R_limit)); */
        force[i] = prefactor * kspring_cov * diff_bond_length[i] / (1.0 - diff_bond_length[i] * \
                                                                    diff_bond_length[i] / (R_limit * R_limit));

        tforce[i] = force[i] * arr_costheta[i];
    }


    /* PRINT HERE */
    /* for ( i = 0; i < m1->num_atoms_ca - 1; i++ ) { */
    /*     printf("ref: %f    curr: %f    diff: %f  pot: %f  force: %f tforce(component): %f \n", \ */
    /*            bond_length_ref[i],bond_length[i],diff_bond_length[i],potential[i],force[i], \ */
    /*            tforce[i]); */
    /* } */


    FILE * fp_tension;
    fp_tension = fopen ("tension_by_residue.dat", "a+");

    // FPRINTF
    fprintf(fp_tension,"# R_limit: %f; kspring_cov %f; prefactor %f;\n",R_limit,kspring_cov,prefactor);
    fprintf(fp_tension,"# from tension analyze: endtoend_magnitude: %f  \n# endtoend_vector: %f %f %f\n", \
            magnitude,endtoend.x,endtoend.y,endtoend.z);

    for ( i = 0; i < m2->num_atoms_ca-1; i++ ) {
        /* Main */
        fprintf(fp_tension,"%6.2f ",tforce[i]);
    }
    fprintf(fp_tension,"\n");
    fclose(fp_tension);
}


void tubulin_monomers_angles( Molecule *m1, Molecule *m2, int r1, int r2, FILE *fp) {
    int i;
    // not yet set
    /* printf("for chain: %d\n",m1->chainid); */
    /* printf("for chain: %d\n",m2->chainid); */

    Vector pos1;
    pos1.x = 0.0;
    pos1.y = 0.0;
    pos1.z = 0.0;
    Vector pos2;
    pos2.x = 0.0;
    pos2.y = 0.0;
    pos2.z = 0.0;

    Vector v1;
    v1.x = 0.0;
    v1.y = 0.0;
    v1.z = 0.0;

    Vector v1_norm;
    v1_norm.x = 0.0;
    v1_norm.y = 0.0;
    v1_norm.z = 0.0;

    Vector pos3;
    pos3.x = 0.0;
    pos3.y = 0.0;
    pos3.z = 0.0;

    Vector pos4;
    pos4.x = 0.0;
    pos4.y = 0.0;
    pos4.z = 0.0;

    Vector v2;
    v2.x = 0.0;
    v2.y = 0.0;
    v2.z = 0.0;

    Vector v2_norm;
    v2_norm.x = 0.0;
    v2_norm.y = 0.0;
    v2_norm.z = 0.0;

    for ( i = 0; i < m1->num_atoms_ca; i++ ) {
        /* printf("%d\n",m1->resid[i]); */
        if ( m1->resid[i] == r1) {
            /* printf(" %f %f %f\n",m1->pos[i].x,m1->pos[i].y,m1->pos[i].z); */
            pos1.x = m1->pos[i].x;
            pos1.y = m1->pos[i].y;
            pos1.z = m1->pos[i].z;
        }
        if ( m1->resid[i] == r2) {
            /* printf(" %f %f %f\n",m1->pos[i].x,m1->pos[i].y,m1->pos[i].z); */
            pos2.x = m1->pos[i].x;
            pos2.y = m1->pos[i].y;
            pos2.z = m1->pos[i].z;
        }
    }
    for ( i = 0; i < m2->num_atoms_ca; i++ ) {
        /* printf("%d\n",m2->resid[i]); */
        if ( m2->resid[i] == r1) {
            /* printf(" %f %f %f\n",m2->pos[i].x,m2->pos[i].y,m2->pos[i].z); */
            pos3.x = m2->pos[i].x;
            pos3.y = m2->pos[i].y;
            pos3.z = m2->pos[i].z;
        }
        if ( m2->resid[i] == r2) {
            /* printf(" %f %f %f\n",m2->pos[i].x,m2->pos[i].y,m2->pos[i].z); */
            pos4.x = m2->pos[i].x;
            pos4.y = m2->pos[i].y;
            pos4.z = m2->pos[i].z;
        }
    }

    get_vector(pos1, pos2, &v1);
    get_vector(pos3, pos4, &v2);
    /* printf("v1: %f %f %f\n",v1.x,v1.y,v1.z); */
    /* printf("v2: %f %f %f\n",v2.x,v2.y,v2.z); */
    normalize(v1, &v1_norm);
    normalize(v2, &v2_norm);
    double costheta;
    double acostheta;
    double deg;
    costheta = get_costheta(v1_norm,v2_norm);
    acostheta = acos(costheta);
    deg = 180.0 * acostheta / M_PI;
    /* printf("angle: %f",costheta); */
    debug("angle: %6.3f  acos: %6.3f  deg: %6.3f\n",costheta,acostheta,deg);
    // FPRINTF
    /* fprintf(fp,"%f ",costheta); */
    fprintf(fp,"%4.2f  ",deg);
}

void get_tubulin_centroid( Molecule *m1, Vector *centroid) {
// void get_tubulin_centroid( Molecule *m1 ) {
// void get_tubulin_centroid( Molecule *m1, std::vector<Vector> *centroid) {
    int i;
    double xtot = 0.0, ytot = 0.0, ztot = 0.0;


    for ( i = 0; i < m1->num_atoms_ca; i++ ) {
        /* printf("%d\n",i); */
        xtot += m1->pos[i].x;
        ytot += m1->pos[i].y;
        ztot += m1->pos[i].z;
    }


    // std::cout << "the_chainid: "<< m1->chainid << std::endl;


    // std::cout << centroid[0] << std::endl;
    // std::cout << centroid[m1->chainid] << std::endl;

    // for(std::vector<Vector>::iterator it = centroid[m1->chainid].begin(); it != centroid[m1->chainid].end(); it++) {
    //     std::cout << it << std::endl;
    // }

    centroid->x = xtot / m1->num_atoms_ca;
    centroid->y = ytot / m1->num_atoms_ca;
    centroid->z = ztot / m1->num_atoms_ca;
}

void get_centroid_from_array( Vector arr_centroid[],int i, int j, Vector *centroid ) {
    double xtot = 0.0, ytot = 0.0, ztot = 0.0;

    debug("centroid from: ||%d %d||\n",i,j);

    for ( int k = i; k <= j; k++ ) {
        /* printf("%f %f %f\n",arr_centroid[k].x,arr_centroid[k].y,arr_centroid[k].z); */
        xtot += arr_centroid[k].x;
        ytot += arr_centroid[k].y;
        ztot += arr_centroid[k].z;
    }

    int t = 0;
    t = j - i + 1;
    // printf("centroid from this many points: %d\n",t);
    centroid->x = xtot / t;
    centroid->y = ytot / t;
    centroid->z = ztot / t;

    debug("centroid_from_array: %f %f %f | %d\n",centroid->x,centroid->y,centroid->z,t);

}

Vector get_centroid_from_centroids(Chain *chain,int start_chain,int stop_chain) {
    double xtot = 0.0, ytot = 0.0, ztot = 0.0;
    Vector centroid;

    debug("centroid from: ||%d %d||\n",start_chain,stop_chain);

    for ( int i=start_chain; i<=stop_chain; i++ ) {
        // printf("%d:%d\n",i,chain[i].chainid);
        xtot += chain[i].centroid.x;
        ytot += chain[i].centroid.y;
        ztot += chain[i].centroid.z;
    }

    int t = 0;
    t = stop_chain - start_chain + 1;
    // printf("centroid from this many points: %d\n",t);
    centroid.x = xtot / t;
    centroid.y = ytot / t;
    centroid.z = ztot / t;

    debug("centroid_of_centroids: %f %f %f | %d\n",centroid.x,centroid.y,centroid.z,t);

    // return chain[start_chain].centroid;
    return centroid;
}
void get_centroid_movement(Chain *chain_ref,Chain *chain_later,int num_chains) {
    int cid;
    // Vector sep1;
    // Vector sep2;
    // Vector movement;
    double d_ref = 0.0;
    double d_later = 0.0;
    double d_move = 0.0;


    FILE * fp_separation_N;
    fp_separation_N = fopen ("centroid_movement_n.dat", "a+");
    FILE * fp_separation_E;
    fp_separation_E = fopen ("centroid_movement_e.dat", "a+");
    FILE * fp_separation_S;
    fp_separation_S = fopen ("centroid_movement_s.dat", "a+");
    FILE * fp_separation_W;
    fp_separation_W = fopen ("centroid_movement_w.dat", "a+");


    for (int i=0; i<num_chains; i++ ) {
        // printf("%d > in contact with > %d %d %d %d\n",i,\
        //        chain_ref[i].LonN,\
        //        chain_ref[i].LatE,\
        //        chain_ref[i].LonS,\
        //        chain_ref[i].LatW);

        for (int j=0; j<4; j++) {

            if (j == 0 ) {
                cid = chain_ref[i].LonN;
            } else if ( j == 1 ) {
                cid = chain_ref[i].LatE;
            } else if ( j == 2 ) {
                cid = chain_ref[i].LonS;
            } else if ( j == 3 ) {
                cid = chain_ref[i].LatW;
            }

            if ( cid == -1 ) {
                // printf("cid: -1 (nil)\n",cid);
                d_ref = 0.0;
                d_later = 0.0;
                d_move = 0.0;
                // continue;
            } else {
                // printf("cid: %d\n",cid);
                d_ref = distance(chain_ref[i].centroid,chain_ref[cid].centroid);
                d_later = distance(chain_later[i].centroid,chain_later[cid].centroid);
                d_move = d_later - d_ref;
            }
            debug("centroid_movement: %f %f %f\n",d_move,d_later,d_ref);

            // printf("%f %f %f\n",chain_ref[i].centroid.x,\
            //        chain_ref[i].centroid.y,\
            //        chain_ref[i].centroid.z);
            // printf("%f %f %f\n",chain_ref[cid].centroid.x,\
            //        chain_ref[cid].centroid.y,\
            //        chain_ref[cid].centroid.z);

            // printf("%f %f %f\n",chain_later[i].centroid.x,\
            //        chain_later[i].centroid.y,\
            //        chain_later[i].centroid.z);
            // printf("%f %f %f\n",chain_later[cid].centroid.x,\
            //        chain_later[cid].centroid.y,\
            //        chain_later[cid].centroid.z);


            if (j == 0 ) {
                fprintf(fp_separation_N," %4.5f",d_move);
            } else if ( j == 1 ) {
                fprintf(fp_separation_E," %4.5f",d_move);
            } else if ( j == 2 ) {
                fprintf(fp_separation_S," %4.5f",d_move);
            } else if ( j == 3 ) {
                fprintf(fp_separation_W," %4.5f",d_move);
            } // if

        } // j 0,1,2,3

    } // i, num_chains

    fprintf(fp_separation_N," \n");
    fprintf(fp_separation_E," \n");
    fprintf(fp_separation_S," \n");
    fprintf(fp_separation_W," \n");

    fclose(fp_separation_N);
    fclose(fp_separation_E);
    fclose(fp_separation_S);
    fclose(fp_separation_W);

    // int LonN; // contacted monomer
    // int LatE;
    // int LonS;
    // int LatW;
    // exit(0);
}
void get_centroid_xyzposition_movement(Chain *chain_later,int used_chains,FILE *fp) {

    for (int i=0; i<used_chains; i++ ) {

        fprintf(fp,"%.2f %.2f %.2f ",chain_later[i].centroid.x, \
                chain_later[i].centroid.y,                      \
                chain_later[i].centroid.z);
    }
    fprintf(fp,"\n");
}

// void print_centroid

void get_latangle_from35centroids(Chain *chain_ref,Chain *chain_later,int m,FILE *fp_ang,int sel=0) {
    // sel: 0: S,N
    // sel: 1: W,E

    //  W:  m1  .. m2 .. [m3] .. m4 .. m5 :E
    // 5
    // int m1,m2,m3,m4,m5;
    // 3
    int m2,m3,m4;

    // --- Neighborhood watch. Is it negative 1? ---
    // debug("W: %d .. %d .. [%d] .. %d .. %d :E\n",m1,m2,m3,m4,m5);
    if (sel == 0) { // S,N
        // m5 = chain_ref[chain_ref[m].LonN].LonN;
        m4 = chain_ref[m].LonN;
        m3 = m;
        m2 = chain_ref[m].LonS;
        // m1 = chain_ref[chain_ref[m].LonS].LonS;

        debug("S: %d .. [%d] .. %d :N\n",m2,m3,m4);
        if ( m2==-1 || m3==-1 || m4==-1) return;

    } else if (sel == 1) { // W,E
        // m5 = chain_ref[chain_ref[m].LatE].LatE;
        m4 = chain_ref[m].LatE;
        m3 = m;
        m2 = chain_ref[m].LatW;
        // m1 = chain_ref[chain_ref[m].LatW].LatW;

        debug("W: %d .. [%d] .. %d :E\n",m2,m3,m4);
        if ( m2==-1 || m3==-1 || m4==-1) return;
    }

    // 5
    // Vector v12,v23,v34,v45;
    // Vector v12n,v23n,v34n,v45n;
    // 3
    Vector v23,v34;
    Vector v23n,v34n;


    //5
    // double cos_vert2,cos_vert3,cos_vert4;
    // cos_vert2 = cos_vert3 = cos_vert4 = 0.0;
    // // debug("cos_vert: 000-> (2)%f (3)%f (4)%f\n",cos_vert2,cos_vert3,cos_vert4);
    // // double cos_vert23, cos_vert34;
    // double vert2,vert3,vert4;
    // vert2 = vert3 = vert4 = 0.0;
    // // debug("vert: 000-> (2)%f (3)%f (4)%f\n",vert2,vert3,vert4);
    // // double vert23,vert34,angle23,angle34;
    // // vert23 = vert34 = 0.0;

    //3
    double cos_vert3;
    cos_vert3 = 0.0;
    // debug("cos_vert: 000-> (2)%f (3)%f (4)%f\n",cos_vert2,cos_vert3,cos_vert4);
    // double cos_vert23, cos_vert34;
    double vert3;
    vert3 = 0.0;
    // debug("vert: 000-> (2)%f (3)%f (4)%f\n",vert2,vert3,vert4);
    // double vert23,vert34,angle23,angle34;
    // vert23 = vert34 = 0.0;


    // 5
    // double angle2,angle3,angle4;
    // angle2 = angle3 = angle4 = 0.0;
    // // debug("angle2-4: 000-> %f %f %f\n",angle2,angle3,angle4);

    // debug("EVAL: W: %d .. %d .. [%d] .. %d .. %d :E\n",m1,m2,m3,m4,m5);
    // // PRINT CENTROIDs
    // debug("centroid1: %f %f %f\n",chain_later[m1].centroid.x,\
    //       chain_later[m1].centroid.y,
    //       chain_later[m1].centroid.z);
    // debug("centroid2: %f %f %f\n",chain_later[m2].centroid.x,\
    //       chain_later[m2].centroid.y,
    //       chain_later[m2].centroid.z);
    // debug("centroid3: %f %f %f\n",chain_later[m3].centroid.x,\
    //       chain_later[m3].centroid.y,
    //       chain_later[m3].centroid.z);
    // debug("centroid4: %f %f %f\n",chain_later[m4].centroid.x,\
    //       chain_later[m4].centroid.y,
    //       chain_later[m4].centroid.z);
    // debug("centroid5: %f %f %f\n",chain_later[m5].centroid.x,\
    //       chain_later[m5].centroid.y,
    //       chain_later[m5].centroid.z);
    // 3
    double angle3;
    angle3 = 0.0;
    // debug("angle2-4: 000-> %f %f %f\n",angle2,angle3,angle4);


    // debug("EVAL: W: %d .. [%d] .. %d :E\n",m2,m3,m4);
    // PRINT CENTROIDs
    // debug("centroid1: %f %f %f\n",chain_later[m1].centroid.x,\
    //       chain_later[m1].centroid.y,
    //       chain_later[m1].centroid.z);
    // debug("centroid2: %f %f %f\n",chain_later[m2].centroid.x,\
    //       chain_later[m2].centroid.y,
    //       chain_later[m2].centroid.z);
    // debug("centroid3: %f %f %f\n",chain_later[m3].centroid.x,\
    //       chain_later[m3].centroid.y,
    //       chain_later[m3].centroid.z);
    // debug("centroid4: %f %f %f\n",chain_later[m4].centroid.x,\
    //       chain_later[m4].centroid.y,
    //       chain_later[m4].centroid.z);
    // debug("centroid5: %f %f %f\n",chain_later[m5].centroid.x,\
    //       chain_later[m5].centroid.y,
    //       chain_later[m5].centroid.z);


    // 5
    // debug("[%d] now getting vectors\n",m);
    // // Vector v12,v23,v34,v45;
    // get_vector(chain_later[m1].centroid,chain_later[m2].centroid,&v12);
    // get_vector(chain_later[m2].centroid,chain_later[m3].centroid,&v23);
    // get_vector(chain_later[m3].centroid,chain_later[m4].centroid,&v34);
    // get_vector(chain_later[m4].centroid,chain_later[m5].centroid,&v45);

    // 3
    // debug("[%d] now getting vectors\n",m);
    // Vector v12,v23,v34,v45;
    // get_vector(chain_later[m2].centroid,chain_later[m3].centroid,&v23);
    v23 = get_vector(chain_later[m2].centroid,chain_later[m3].centroid);
    // get_vector(chain_later[m3].centroid,chain_later[m4].centroid,&v34);
    v34 = get_vector(chain_later[m3].centroid,chain_later[m4].centroid);

    // ONLY FOR EAST/WEST, to put latitudinal monomers in a plane.
    // if (sel==1) {
    //     // zeroing z component
    //     // 5
    //     // v12.z = v23.z = v34.z = v45.z = 0.0;
    //     // 3
    //     v23.z = v34.z = 0.0;
    // }


    // 5
    // // double normalize ( Vector v1, Vector *v2 );
    // // Vector v12n,v23n,v34n,v45n;
    // normalize(v12,&v12n);
    // normalize(v23,&v23n);
    // normalize(v34,&v34n);
    // normalize(v45,&v45n);
    // 3
    normalize(v23,&v23n);
    normalize(v34,&v34n);


    // PRINT NORMALIZED
    // printf("angles of interest: v12n - v34n;   v23n - v45n\n"); // NO
    // printf("v12n: %f %f %f\n",v12n.x,v12n.y,v12n.z);
    // printf("v23n: %f %f %f\n",v23n.x,v23n.y,v23n.z);
    // printf("v34n: %f %f %f\n",v34n.x,v34n.y,v34n.z);
    // printf("v45n: %f %f %f\n",v45n.x,v45n.y,v45n.z);


    // 5
    // // double cos_vert2,cos_vert3,cos_vert4;
    // cos_vert2 = get_costheta(v12n,v23n);
    // cos_vert3 = get_costheta(v23n,v34n);
    // cos_vert4 = get_costheta(v34n,v45n);
    // 3
    // double cos_vert2,cos_vert3,cos_vert4;
    cos_vert3 = get_costheta(v23n,v34n);


    // cos_vert23 = get_costheta(v12n,v34n);
    // cos_vert34 = get_costheta(v23n,v45n);

    // double vert2,vert3,vert4;
    // if ( )
    // 5
    // vert2 = acos(cos_vert2);
    // vert3 = acos(cos_vert3);
    // vert4 = acos(cos_vert4);
    // 3
    vert3 = acos(cos_vert3);

    // double vert23,vert34;
    // vert23 = acos(cos_vert23);
    // vert34 = acos(cos_vert34);

    // 5
    // angle2 = vert2 * 180.0 / M_PI;
    // angle3 = vert3 * 180.0 / M_PI;
    // angle4 = vert4 * 180.0 / M_PI;
    // 3
    angle3 = vert3 * 180.0 / M_PI;


    // angle23 = vert23 * 180.0 / 3.1415;
    // angle34 = vert34 * 180.0 / 3.1415;

    // printf("monomer_angles: (2):%f  (3):%f  (4):%f\n",cos_vert2 * 180 / 3.14,
    //        cos_vert3 * 180.0 / 3.14,
    //        cos_vert4 * 180.0 / 3.14);
    // printf("monomer_angles (cos):\t (2):%7.4f  (3):%7.4f  (4):%7.4f\n",cos_vert2,cos_vert3,cos_vert4);
    // printf("monomer_angles (acos):\t (2):%7.4f  (3):%7.4f  (4):%7.4f\n",vert2,vert3,vert4);
    // printf("monomer_angles (deg):\t (2):%7.4f  (3):%7.4f  (4):%7.4f\n",vert2 * 180 / 3.14,vert3 * 180 / 3.14,vert
    // 4 * 180 / 3.14);

    // printf("intermonomer_angles: (cos):\t (23)%7.4f  (34)%7.4f\n",cos_vert23,cos_vert34);
    // printf("intermonomer_angles: (acos):\t (23)%7.4f  (34)%7.4f\n",vert23,vert34);
    // printf("intermonomer_angles: (deg):\t (23)%7.4f  (34)%7.4f\n",angle23,angle34);
// }
// printf("m: %d || west: %f  center: %f  east: %f\n",m,vert2,vert3,vert4);


// printf("\n");
// monomer West East Angle   monomer2 W2 E2 Angle2 monomer4 W4 E4 Angle4
// m3      m2   m4   angle3  m2       m1 m3 angle2 m4       m3 m5 angle4


    // 5
    // if (sel==0) {
    //     FILE * fp_angles_ns;
    //     fp_angles_ns = fopen ("mt_angles_ns.dat", "a+");
    //     fprintf(fp_angles_ns, "%4d %4d %4d  %7.3f    ",m3,m2,m4,angle3);
    //     fprintf(fp_angles_ns, "%4d %4d %4d  %7.3f    ",m2,m1,m3,angle2);
    //     fprintf(fp_angles_ns, "%4d %4d %4d  %7.3f\n",m4,m3,m5,angle4);
    //     fclose(fp_angles_ns);
    // } else if (sel==1){
    //     FILE * fp_angles_ew;
    //     fp_angles_ew = fopen ("mt_angles_ew.dat", "a+");
    //     fprintf(fp_angles_ew, "%4d %4d %4d  %7.3f    ",m3,m2,m4,angle3);
    //     fprintf(fp_angles_ew, "%4d %4d %4d  %7.3f    ",m2,m1,m3,angle2);
    //     fprintf(fp_angles_ew, "%4d %4d %4d  %7.3f\n",m4,m3,m5,angle4);
    //     fclose(fp_angles_ew);
    // }
    // 3



    // if (sel==0) {
    //     FILE * fp_angles_ns;
    //     fp_angles_ns = fopen ("mt_angles_ns.dat", "a+");
    //     fprintf(fp_angles_ns, "%4d %4d %4d  %7.3f\n",m3,m2,m4,angle3);
    //     // fprintf(fp_angles_ns, "%4d %4d %4d  %7.3f    ",m2,m1,m3,angle2);
    //     // fprintf(fp_angles_ns, "%4d %4d %4d  %7.3f\n",m4,m3,m5,angle4);
    //     fclose(fp_angles_ns);
    // } else if (sel==1){
    //     FILE * fp_angles_ew;
    //     fp_angles_ew = fopen ("mt_angles_ew.dat", "a+");
    //     fprintf(fp_angles_ew, "%4d %4d %4d  %7.3f\n",m3,m2,m4,angle3);
    //     // fprintf(fp_angles_ew, "%4d %4d %4d  %7.3f    ",m2,m1,m3,angle2);
    //     // fprintf(fp_angles_ew, "%4d %4d %4d  %7.3f\n",m4,m3,m5,angle4);
    //     fclose(fp_angles_ew);
    // }

    fprintf(fp_ang,"%3d %3d %3d %6.2f   ",m3,m2,m4,angle3);

    return;
}

// std::vector<double> get_3angles_dimerbyalpha(Chain *chain_ref,Chain *chain_later,int m,std::vector<double> lst_angles3)
std::vector<double> get_3angles_dimerbyalpha(Chain *chain_ref,Chain *chain_later,int cid)
{
    int alpha,beta,alpha_w,alpha_e,beta_w,beta_e,alpha2;

    // 3 angles: int m is A1. the alpha monomer.
    //    ^ w      w ^
    //   1. A1 --  B .2 -- A2:   < - 3 - >
    //    v e      e v

    alpha = cid;
    beta = chain_ref[alpha].LonN;
    alpha_w = chain_ref[alpha].LatW;
    alpha_e = chain_ref[alpha].LatE;
    beta_w = chain_ref[beta].LatW;
    beta_e = chain_ref[beta].LatE;
    alpha2 = chain_ref[beta].LonN;


    // Centroids:
    chain_later[alpha].ComputeCentroid();
    chain_later[beta].ComputeCentroid();

    if (alpha_w != -1)
    {
        chain_later[alpha_w].ComputeCentroid();
    }
    if (alpha_e != -1)
    {
        chain_later[alpha_e].ComputeCentroid();
    }
    if (beta_w != -1)
    {
        chain_later[beta_w].ComputeCentroid();
    }
    if (beta_e != -1)
    {
        chain_later[beta_e].ComputeCentroid();
    }
    if (alpha2 != -1)
    {
        chain_later[alpha2].ComputeCentroid();
    }

    // Vector
    // one: WalphaE
    // two: WbetaE
    // three: alphabetaN
    // ba1: a1 - b (backward), a2 - b (forward);
    Vector one_w,one_e,two_w,two_e,three_ba1,three_ba2;
    double vert1,vert2,vert3,cos_vert1,cos_vert2,cos_vert3;
    double angle1,angle2,angle3;
    vert1 = vert2 = vert3 = 0.0;
    cos_vert1 = cos_vert2 = cos_vert3 = 0.0;
    angle1 = angle2 = angle3 = 0.0;


    // vectors:
    if ((alpha_w != -1) && (alpha_e != -1))
    {
        one_w = get_vector(chain_later[alpha].centroid,chain_later[alpha_w].centroid);
        one_e = get_vector(chain_later[alpha].centroid,chain_later[alpha_e].centroid);
    }
    if ((beta_w != -1) && (beta_e != -1))
    {
        two_w = get_vector(chain_later[beta].centroid,chain_later[beta_w].centroid);
        two_e = get_vector(chain_later[beta].centroid,chain_later[beta_e].centroid);
    }
    if (alpha2 != -1)
    {
        three_ba1 = get_vector(chain_later[beta].centroid,chain_later[alpha].centroid);
        three_ba2 = get_vector(chain_later[beta].centroid,chain_later[alpha2].centroid);
    }


    // normalize && get_costheta && get cos(radians) && degrees.
    Vector norm_one_w,norm_one_e,norm_two_w,norm_two_e,norm_three_ba1,norm_three_ba2;

    if ((alpha_w != -1) && (alpha_e != -1))
    {
        norm_one_w = normalize(one_w);
        norm_one_e = normalize(one_e);
        cos_vert1 = get_costheta(norm_one_w,norm_one_e);
        vert1 = acos(cos_vert1);
        angle1 = vert1 * 180.0 / M_PI;
    }
    if ((beta_w != -1) && (beta_e != -1))
    {
        norm_two_w = normalize(two_w);
        norm_two_e = normalize(two_e);
        cos_vert2 = get_costheta(norm_two_w,norm_two_e);
        vert2 = acos(cos_vert2);
        angle2 = vert2 * 180.0 / M_PI;
    }
    if (alpha2 != -1)
    {
        norm_three_ba1 = normalize(three_ba1);
        norm_three_ba2 = normalize(three_ba2);
        cos_vert3 = get_costheta(norm_three_ba1,norm_three_ba2);
        vert3 = acos(cos_vert3);
        angle3 = vert3 * 180.0 / M_PI;
    }

    std::vector<double> lst_angles3;
    lst_angles3.push_back(angle1);
    lst_angles3.push_back(angle2);
    lst_angles3.push_back(angle3);

    return lst_angles3;
} // get_3angles_dimerbyalpha


// void get_midpoint_reference_indices(Chain *chain_ref,int *index_i,int *index_c,int sel) {
// void get_midpoint_reference_indices(Chain *chain_ref,int i,int *index_i,int *index_c,int sel) {

//     int cid;
//     // double d_ref = 0.0;
//     Vector vec_centroids;
//     Vector v_midpoint;

//     if (sel == 1 ) {
//         cid = chain_ref[i].LonN;
//     } else if ( sel == 2 ) {
//         cid = chain_ref[i].LatE;
//     } else if ( sel == 3 ) {
//         cid = chain_ref[i].LonS;
//     } else if ( sel == 4 ) {
//         cid = chain_ref[i].LatW;
//     }

//     // d_ref = distance(chain_ref[i].centroid,chain_ref[cid].centroid);
//     get_vector(chain_ref[i].centroid,chain_ref[cid].centroid,&vec_centroids);
//     debug("centroid_vector: %f %f %f",vec_centroids.x,vec_centroids.y,vec_centroids.z);
//     // midpoint of vector between centroids
//     v_midpoint = midpoint(vec_centroids);
//     debug("centroid_midpoint: %f %f %f",v_midpoint.x,v_midpoint.y,v_midpoint.z);
//     // get index closest to point


//     exit(0);
// }

void get_monomer_monomer_info( Molecule *m1, Molecule *m2,\
                               Vector ref_norm, Vector m1centroid,\
                               Vector m2centroid ) {

    /* printf(" %f %f %f\n",ref_norm.x,ref_norm.y,ref_norm.z); */
    /* printf(" %d in contact with %d \n",m1->num_atoms_ca,m2->num_atoms_ca); */
    /* printf(" %f %f %f -- %f %f %f\n",m1centroid.x,m1centroid.y,m1centroid.z,\ */
    /*        m2centroid.x,m2centroid.y,m2centroid.z); */

    Vector param_v12;
    get_vector(m1centroid, m2centroid, &param_v12);
    Vector param_v12n;
    normalize(param_v12, &param_v12n);
    /* printf(" %f %f %f\n",param_v12n.x,param_v12n.y,param_v12n.z); */


    /* CHECK */
    /* Vector centroid_check1; */
    /* Vector centroid_check2; */
    /* Vector v12; */
    /* get_tubulin_centroid(m1, &centroid_check1); */
    /* get_tubulin_centroid(m2, &centroid_check2); */
    /* get_vector(centroid_check1, centroid_check2, &v12); */
    /* Vector normcheck; */
    /* /\* printf(" %f %f %f\n",centroid_check1.x,centroid_check1.y,centroid_check1.z); *\/ */
    /* /\* printf(" %f %f %f\n",centroid_check1n.x,centroid_check1n.y,centroid_check1n.z); *\/ */
    /* normalize(v12, &normcheck); */
    /* printf(" %f %f %f\n",normcheck.x,normcheck.y,normcheck.z); */

}

// void compute_nearest_neighbors(Chain *chain, int num_chains) {
void compute_nearest_neighbors(Chain *arr, int num_chains) {

    debug("Inside compute_nearest_neigbhors function(chain.cpp); using 132.0 neighbor cutoff.\n");

    double dist;
    int count;

    for ( int i=0; i < num_chains; i++ ) {
        // printf("ichain: %d\n",i);
        count = 0;
        for ( int j=0; j < num_chains; j++ ) {
            // printf("%d %d\n",j,arr[j].chainid);
            if ( i == j ) {
                continue;
            } else {
                dist = distance(arr[i].centroid,arr[j].centroid);
                if ( dist < 132.0 ) {
                    //get broken at != -1
                    for ( int k=0; k < 32; k++ ) {
                        if ( arr[i].nearest_targets[k] == -1 ) {
                            arr[i].nearest_targets[k] = j;
                            arr[i].nearest_targets_dist[k] = dist;
                            break;
                        }
                    }
                    count++;

                    // PRINT HERE (targets)!
                    // printf("%d-%d  (c:%d d:%f)\n",i,j,count,dist);
                }
            }

        }
    }

    // PRINT HERE: TARGETS
    // for ( int i=0; i < num_chains; i++ ) {
    //     for ( int k=0; k < 32; k++ ) {
    //         if ( arr[i].nearest_targets_dist[k] != 0.0 ){
    //             // printf(" %d %f\n",k,arr[i].nearest_targets_dist[k]);
    //         }
    //     }
    // }

    double min_dist;
    int index;

    // get the lowest (first)
    for ( int i=0; i < num_chains; i++ ) {
        min_dist = 133.0;
        for ( int k=0; k < 32; k++ ) {
            if ((arr[i].nearest_targets_dist[k] < min_dist) && \
                (arr[i].nearest_targets_dist[k] != 0.0)) {
                index = arr[i].nearest_targets[k];
                min_dist = arr[i].nearest_targets_dist[k];
                // printf("%d %f\n",index,min_dist);
            }
        }
        arr[i].nearestneighbors_chainid[0] = index;
        arr[i].nearestneighbors_dist[0] = min_dist;
    }

    // PRINT HERE: lowest or position in the nearestneighbors_dist array [size 8]
    // for ( int i=0; i < num_chains; i++ ) {
    //     printf("%d %f\n",i,arr[i].nearestneighbors_dist[0]);
    // }
    // exit(0);


    for ( int i=0; i < num_chains; i++ ) {
        for ( int e=1; e < 8; e++) {
            min_dist = 133.0;
            for ( int k=0; k < 32; k++ ) {
                if (arr[i].nearest_targets_dist[k] < min_dist && \
                    arr[i].nearest_targets_dist[k] > arr[i].nearestneighbors_dist[e-1]) {
                    index = arr[i].nearest_targets[k];
                    min_dist = arr[i].nearest_targets_dist[k];
                    // printf("%d %f\n",index,min_dist);
                }
            }
            arr[i].nearestneighbors_chainid[e] = index;
            arr[i].nearestneighbors_dist[e] = min_dist;
        }
    }

    // PRINT HERE: results ..
    // for ( int i=0+59; i < num_chains; i++ ) {
    // for ( int i=0+5; i < 10; i++ ) {
    //     printf("%d ",i);
    //     printf("indices: %d|%d\n",arr[i].index,arr[i].findex);
    //     for ( int e=0; e < 8; e++) {
    //         printf(" %d|%.5f", \
    //                arr[i].nearestneighbors_chainid[e],\
    //                arr[i].nearestneighbors_dist[e]);
    //     }
    //     printf("\n");
    //     for ( int e=0; e < 8; e++) {
    //         printf(" [%d %d]", \
    //                arr[arr[i].nearestneighbors_chainid[e]].index,\
    //                arr[arr[i].nearestneighbors_chainid[e]].findex);
    //     }
    //     // printf("\nindices: %d %d",arr[i].index,arr[i].findex);
    //     printf("\n");
    // }
    // exit(0);



    // double min_dist;
    // double keep_dist = 133.0;
    // int keep_k;
    // int e_value;
    // int pos_nearest;

    // // Evaluate the targets.
    // for ( int i=0; i < num_chains; i++ ) {
    //     // min_dist = 0.0;
    //     // keep_dist = 133.0;
    //     printf("CHAIN %d\n",i);
    //     for ( int ei=0; ei < 8; ei++ ) {
    //         // printf("ei: %d\n");
    //     // while ( e_value < 8 ) {
    //         for ( int k=0; k < 32; k++ ) {
    //             if ( arr[i].nearest_targets[k] != -1 ) {
    //                 printf("%d %d\n",k,arr[i].nearest_targets[k]);
    //                 // continue;
    //             } else {
    //                 break;
    //             }

    //             for ( int e=0; e < 8; e++ ) {
    //                 printf("ei: %d   E: %d  MIN:%f\n",ei,e,min_dist);
    //                 dist = distance(arr[i].centroid,arr[arr[i].nearest_targets[e]].centroid);
    //                 printf("comparing: %d with %d, distaway: %f\n",arr[i].nearest_targets[e],i,dist);

    //                 if (( arr[i].nearestneighbors_dist[0] == 0.0) && ( e==0 )) {
    //                     printf("1st case\n");
    //                     // min_dist = 0.0;
    //                     e_value = e;
    //                     break;
    //                 } else if ((arr[i].nearestneighbors_dist[e] == 0.0) && (e != 0)) {
    //                     printf("case 2-7\n");
    //                     pos_nearest = e - 1;
    //                     min_dist = arr[i].nearestneighbors_dist[pos_nearest];
    //                     printf("assigning MIN_dist: %f pos_nearest: %d\n",min_dist,pos_nearest);
    //                     e_value = e;
    //                     break;
    //                 } else {
    //                     printf("3rd case\n");
    //                     e_value = 8;
    //                     // exit(0);
    //                     break;
    //                 }


    //             //     //     if ( arr[i].nearestneighbors_dist[e] == 0.0 ) {
    //             //     //     pos_nearest = e - 1;
    //             //     //     min_dist = arr[i].nearestneighbors_dist[pos_nearest];
    //             //     //     printf("assigning min_dist: %f",min_dist);
    //             //     //     e_value = e;
    //             //     // } else {
    //             //     //     printf("say what?\n");
    //             //     //     exit(1);
    //             //     // }
    //             // } // entry e 0-7 (8)
    //             // printf("e_value: %d   min_dist: %f\n",e_value,min_dist);
    //             // dist = distance(arr[i].centroid,arr[arr[i].nearest_targets[k]].centroid);
    //             // printf("current dist: %f  vs.  %f <- keep_dist\n",dist,keep_dist);
    //             // if (( dist > min_dist ) && ( dist < keep_dist) ){
    //             //     printf("changing keep_dist... %f\n",dist);
    //             //     keep_dist = dist;
    //             //     keep_k = k;
    //             //     // break; //
    //             }
    //             printf("track: ei:%d k:%d pos_nearest:%d e_value:%d min_dist:%f");
    //             double min_dist;
    //             double keep_dist = 133.0;
    //             int keep_k;
    //             int e_value;
    //             int pos_nearest;

    //         } // targets k 0 - 31 (32)

    //         // exit(0);

    //         // printf("  found %d %f %d\n",keep_k,keep_dist,e_value);
    //         // arr[i].nearestneighbors_chainid[e_value] = keep_k;
    //         // arr[i].nearestneighbors_dist[e_value] = keep_dist;
    //         // break;

    //         } // first for e 8
    //         exit(0);

    //     // } // while
    // } // centroids i


    // check num_chains
    // for ( int i=0; i < 1; i++ ) {
    //     printf("%d\n",i);
    //     for ( int e=0; e < 8; e++ ) {
    //         printf("nearestneighbors: %d %f    ",arr[i].nearestneighbors_chainid[e],\
    //                arr[i].nearestneighbors_dist[e]);
    //     }
    //     printf("\n");
    // }


    //     if ( arr->chainid == j ) {
    //         // printf("\n match. \n");
    //         continue;
    //     } else {
    //         // printf("%d\n",j);
    //         // dist = distance(chain->centroid,)
    //     }
    // }
}

void determine_latlon_neighbors(Chain *chains, int num_chains, Vector axis) {

    double costheta,sintheta;
    Vector cc; // centroid to centroid;
    Vector cc_norm; // normalize centroid to centroid vector;
    double dist;

    // // 20%   0...20...*...80...100.
    // double min_z;
    // double max_z;
    // printf("finding the minimum and maximum on z ..\n");
    // for ( int i=0; i < num_chains; i ++ ) {
    //     if ( i == 0 ) {
    //         min_z = chains[i].centroid.z;
    //         max_z = chains[i].centroid.z;
    //     } else {
    //         if ( min_z > chains[i].centroid.z ) {
    //             min_z = chains[i].centroid.z;
    //         }
    //         if ( max_z < chains[i].centroid.z ) {
    //             max_z = chains[i].centroid.z;
    //         }
    //     }
    // }
    // printf("the minimum z(%f) and maximum z(%f) have been found\n",min_z,max_z);
    // double z_thres = (max_z - min_z) * 0.30;
    // double min_z_shift = min_z + z_thres;
    // double max_z_shift = max_z - z_thres;
    // printf("the threshold minimum z(%f) and maximum z(%f) have been determined\n",min_z_shift,max_z_shift);


    for ( int i=0; i < num_chains; i ++ ) {
        // printf("%d\n",chains[i].chainid);
        // printf(" %f %f %f\n",chains[i].centroid.x,chains[i].centroid.y,chains[i].centroid.z);

        // if (( chains[i].centroid.z < min_z_shift) || ( chains[i].centroid.z > max_z_shift )) {
        //         continue;
        // }

        for ( int j=0; j < 4; j++ ) {
            // printf(" %d ",chains[i].nearestneighbors_chainid[j]);

            get_vector(chains[i].centroid,\
                       chains[chains[i].nearestneighbors_chainid[j]].centroid,\
                       &cc);
            normalize(cc, &cc_norm);
            costheta = get_costheta(axis,cc_norm);
            sintheta = get_sintheta(axis,cc_norm);
            dist = distance(chains[i].centroid,chains[chains[i].nearestneighbors_chainid[j]].centroid);

            // PRINT HERE
            // printf("[%f|%f]",costheta,sintheta);
            // printf(" %f ",dist);
            // printf(">> %f %f %f ",chains[chains[i].nearestneighbors_chainid[j]].centroid.x,\
            //        chains[chains[i].nearestneighbors_chainid[j]].centroid.y, \
            //        chains[chains[i].nearestneighbors_chainid[j]].centroid.z);
            // printf("\n");

            // evaluate the latitudinal / longitudinal comparison.
            if ( costheta > 0.95 ) {
                chains[i].LonN = chains[i].nearestneighbors_chainid[j];
            } else if ( costheta < -0.95 ) {
                chains[i].LonS = chains[i].nearestneighbors_chainid[j];
            }

            if ( sintheta > 0.94 ) {
                if ( costheta > 0.0 ) {
                    chains[i].LatW = chains[i].nearestneighbors_chainid[j];
                } else if ( costheta < 0.0 ) {
                    chains[i].LatE = chains[i].nearestneighbors_chainid[j];
                } else if ( costheta == 0.0 ) {
                    printf("warning: costheta hit 0.0\n");
                    exit(1);
                }
            }
        }

        // PRINT HERE
        // printf("N:%d S:%d W:%d E:%d\n",chains[i].LonN,chains[i].LonS,chains[i].LatW,chains[i].LatE);
        // printf("\n");
    }
}


// Contact build_contact_map_inter(Chain *chain,Contact *map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS], int sel) {
// Contact build_contact_map_inter(Chain *chain,Contact (*map)[MAX_CONTACTS], int sel) {
// Contact build_contact_map_inter(Chain *chain,Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS], int sel) {
// void build_contact_map_inter(Chain *chain,int cid, Contact *map[][MAX_CONTACTS], int sel) {
void build_contact_map_inter_monomer(Chain *chain,int cid, Contact (*map)[MAX_CONTACTS], int sel) {
    // *chain (all of the chains)
    // cid is the selected chain.
    // sel: N E S W | 1 2 3 4
    int count = 0;
    int total = 0;
    int inter_chain;
    double dist;

    switch(sel) {
    case 1:
        inter_chain = chain[cid].LonN;
        break;
    case 2:
        inter_chain = chain[cid].LatE;
        break;
    case 3:
        inter_chain = chain[cid].LonS;
        break;
    case 4:
        inter_chain = chain[cid].LatW;
        break;
    default:
        std::cout << "Invalid Selection. Please try Again." << std::endl;
    } //switch
    if ( inter_chain == -1 ) {
        return;
    }
    // debug("[%d] with %d %d %d %d",cid,chain[cid].LonN,chain[cid].LatE,chain[cid].LonS,chain[cid].LatW);
    // debug("[%d] with <%d> selected",cid,inter_chain);


#ifdef INFO
    printf("evaluating %d__%d <->  %d  NESW(%d)\n",chain[cid].chainid,cid,inter_chain,sel);
#endif // INFO

    // if ( sel==2 ) {

    // }

    for( int i=0; i<chain[cid].num_atoms_ca; i++ ) {
        // for( int j=0; j<chain[inter_chain].num_atoms_ca; j++) {



        for( int j=0; j<chain[inter_chain].num_atoms_ca; j++) {
        // for( int j=i; j<chain[inter_chain].num_atoms_ca; j++) {
            // printf("the count was: %d\n", count);

            if ( i >= j-2 && i <= j+2) {
                // exclude j = i +/- 2
                /* printf("%d\t no --> %d j: %d %d\n",i,j,j-2,j+2); */
            } else {
                /* printf("%d\n",j); */
                dist = distance( chain[cid].pos[i], chain[inter_chain].pos[j]);
                // printf("%8.4f\n",dist);

                if ( dist <= 8.0 ) {
                // if ( dist <= 20.0 ) {
                    /* printf("criterion-> %d %d %8.4f\n",i,j,dist); */


                    /* SINGLE */
                    /* for( count = 0; count < 40; count++ ) { */
                    /*     if ( mol->contact_map[i][count].cresid != -1 ) { */
                    /*         continue; */
                    /*     } else { */
                    /*         mol->contact_map[i][count].cresid = j; */
                    /*         mol->contact_map[i][count].distance = dist; */
                    /*         break; */
                    /*     } */
                    /* } */

                    /* DOUBLE */
                    for( int k=0; k<MAX_CONTACTS; k++ ) {

                        if ( map[i][k].cresid != -1 ) {
                            continue;
                        } else {
                            map[i][k].cresid = j;
                            map[i][k].distance = dist;
                            count ++;
                            break;
                        }
                    }

                    /* and j */
                    for( int l=0; l<MAX_CONTACTS; l++ ) {

                        if (map[j][l].cresid != -1 ) {
                            continue;
                        } else {
                            map[j][l].cresid = i;
                            map[j][l].distance = dist;
                            count ++;
                            break;
                        }
                    }
                    // count ++;

                } /* END 8.0 */
            } // else
        } // j
        // printf("count_curr: %d in %d\n",count,i);
        total += count;
        count = 0;
    } // i
    // printf("the count was: %d\n", count);

    // PRINT HERE
    // debug("inter-total NESW(%d): ch %d, neighbor %d,  tot:%d\n",sel,cid,inter_chain,total);


    switch(sel) {
    case 1:
        chain[cid].LonN_total2x = total;
        break;
    case 2:
        chain[cid].LatE_total2x = total;
        break;
    case 3:
        chain[cid].LonS_total2x = total;
        break;
    case 4:
        chain[cid].LatW_total2x = total;
        break;
    default:
        std::cout << "Invalid Selection. Please try Again." << std::endl;
    } //switch
    // if ( inter_chain == -1 ) {
    //     return;
    // }


    return;
}

void build_contact_map_inter(Chain *chain,int cid1,int cid2,Contact (*map)[MAX_CONTACTS]) {
    // *chain (all of the chains)
    // cid1 is the selected chain.
    // in contact with cid2.
    int count = 0;
    int total = 0;
    double dist;


    if (cid2 == -1) {
        debug("returning immediately(cid) %d <-> %d <---because negative one\n",cid1,cid2);
        return;
    }

    // PRINT HERE
    // debug("evaluating(id)  %d[ca:%d] <-> %d[ca:%d]",chain[cid1].chainid, \
    //       chain[cid1].num_atoms_ca,\
    //       chain[cid2].chainid,\
    //       chain[cid2].num_atoms_ca);


    for( int i=0; i<chain[cid1].num_atoms_ca; i++ ) {
        // for( int j=0; j<chain[cid2].num_atoms_ca; j++) {

        // for( int j=i; j<chain[cid2].num_atoms_ca; j++) {
        for( int j=0; j<chain[cid2].num_atoms_ca; j++) {
            // printf("the count was: %d\n", count);

            if ( i >= j-2 && i <= j+2) {
                // exclude j = i +/- 2
                /* printf("%d\t no --> %d j: %d %d\n",i,j,j-2,j+2); */
            } else {
                // printf("%d\n",j);
                dist = distance( chain[cid1].pos[i], chain[cid2].pos[j]);
                // printf("%8.4f\n",dist);
                // debug("%8.4f",dist);

                // UNDEBUG HERE FIRST
                debug("c:%d-%d\n",i,j);

                if ( dist <= 8.0 ) {
                    /* printf("criterion-> %d %d %8.4f\n",i,j,dist); */

                    // PRINT HERE
                    debug("FOUND ONE: %8.4f __ %d-%d\n",dist,i,j);

                    /* SINGLE */
                    /* for( count = 0; count < 40; count++ ) { */
                    /*     if ( mol->contact_map[i][count].cresid != -1 ) { */
                    /*         continue; */
                    /*     } else { */
                    /*         mol->contact_map[i][count].cresid = j; */
                    /*         mol->contact_map[i][count].distance = dist; */
                    /*         break; */
                    /*     } */
                    /* } */

                    /* DOUBLE */
                    for( int k=0; k<MAX_CONTACTS; k++ ) {

                        if ( map[i][k].cresid != -1 ) {
                            continue;
                        } else {
                            map[i][k].cresid = j;
                            map[i][k].distance = dist;
                            map[i][k].index = i + chain[cid1].index;
                            map[i][k].cindex = j + chain[cid2].index;
                            count ++;
                            break;
                        }
                    }

                    /* and j */
                    for( int l=0; l<MAX_CONTACTS; l++ ) {

                        if (map[j][l].cresid != -1 ) {
                            continue;
                        } else {
                            map[j][l].cresid = i;
                            map[j][l].distance = dist;
                            map[j][l].index = j + chain[cid2].index;
                            map[j][l].cindex = i + chain[cid1].index;
                            count ++;
                            break;
                        }
                    }
                    // count ++;

                } /* END 8.0 */
            } // else
        } // j
        // printf("count_curr: %d in %d\n",count,i);
        total += count;
        // count = 0;
    } // i

    // PRINT HERE
    debug("the count was: %d\n", count);
    debug("inter-total(%d): tot:%d\n",cid1,total);

    return;
}

void compare_contacts(Chain *chref,Chain *chain,int cid, int sel) {
    // and write get chain[i].LatE_total2x;

    int inter_chain = 0;
    Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    Contact refmap[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    int match_count = 0;

    // printf("map >> ");
    // print_contact_map(map);
    // printf("refmap >> ");
    // print_contact_map(refmap);


    switch(sel) {
    case 1:
        inter_chain = chref[cid].LonN;
        memcpy(&map,&chain[cid].contactsLonN, sizeof(chain[cid].contactsLonN));
        // printf("map 1>> \n");
        // print_contact_map(chain[cid].contactsLonN);
        // print_contact_map(map);
        memcpy(&refmap,&chref[cid].contactsLonN, sizeof(chref[cid].contactsLonN));
        // printf("refmap 1>> \n");
        // print_contact_map(refmap);
        break;
    case 2:
        inter_chain = chref[cid].LatE;
        memcpy(&map,&chain[cid].contactsLatE, sizeof(chain[cid].contactsLatE));
        memcpy(&refmap,&chref[cid].contactsLatE, sizeof(chref[cid].contactsLatE));
        break;
    case 3:
        inter_chain = chref[cid].LonS;
        memcpy(&map,&chain[cid].contactsLonS, sizeof(chain[cid].contactsLonS));
        // printf("map 1>> \n");
        // print_contact_map(chain[cid].contactsLonS);
        // print_contact_map(map);
        memcpy(&refmap,&chref[cid].contactsLonS, sizeof(chref[cid].contactsLonS));
        // printf("refmap 1>> \n");
        // print_contact_map(chref[cid].contactsLonS);
        // print_contact_map(refmap);
        // memcpy(map, chain[cid].contactsLonS, sizeof(chain[cid].contactsLonS));
        // memcpy(refmap, chref[cid].contactsLonS, sizeof(chref[cid].contactsLonS));
        break;
    case 4:
        inter_chain = chref[cid].LatW;
        memcpy(&map,&chain[cid].contactsLatW, sizeof(chain[cid].contactsLatW));
        memcpy(&refmap,&chref[cid].contactsLatW, sizeof(chref[cid].contactsLatW));
        break;
    default:
        std::cout << "Invalid Selection. Please try Again." << std::endl;
    } //switch
    if ( inter_chain == -1 ) {
        // exit(1);
        debug("no interface found! (%d)\n",sel);
        return;
    }

    // printf("printing maps .....................\n");
    // print_contact_map(map);
    // print_contact_map(refmap);
    for ( int i=0; i<MOLECULE_INITIAL_SIZE; i++ ) {
        for ( int j=0; j<MAX_CONTACTS; j++ ) {
            if ( map[i][j].cresid == -1 ) {
                continue;
            }
            // printf("%d %f\n",map[i][j].cresid,map[i][j].distance);
            // int atom = i;
            int first = i;
            int second = map[i][j].cresid;

            for ( int k=0; k<MAX_CONTACTS; k++ ) {
                if ( refmap[first][k].cresid == second ) {
                    match_count++;
                    break;
                }
                if ( refmap[second][k].cresid == first ) {
                    match_count++;
                    break;
                }
            } // for3
        } // for2
    } // for1

    // PRINT HERE
    // debug("printing reference map\n");
    // print_contact_map(refmap);
    // debug("printing contact map\n");
    // print_contact_map(map);
    // debug("match_total: %d\n",match_count);


    switch(sel) {
    case 1:
        chain[cid].LonN_total2x = match_count;
        break;
    case 2:
        chain[cid].LatE_total2x = match_count;
        break;
    case 3:
        chain[cid].LonS_total2x = match_count;
        break;
    case 4:
        chain[cid].LatW_total2x = match_count;
        break;
    default:
        std::cout << "Invalid Selection. Please try Again." << std::endl;
    } //switch
}

void evaluate_original_contacts_now(Chain *chref,Chain *chain,int cid, int sel) {
    // return;
    // and write get chain[i].LatE_total2x;

    int inter_chain = 0;
    // Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    Contact refmap[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    int match_count = 0;
    // print_contact_map(refmap);


    switch(sel) {
    case 1:
        inter_chain = chref[cid].LonN;
        // memcpy(&map,&chain[cid].contactsLonN, sizeof(chain[cid].contactsLonN));
        // printf("map 1>> \n");
        // print_contact_map(chain[cid].contactsLonN);
        // print_contact_map(map);
        memcpy(&refmap,&chref[cid].contactsLonN, sizeof(chref[cid].contactsLonN));
        // printf("refmap 1>> \n");
        // print_contact_map(refmap);
        break;
    case 2:
        inter_chain = chref[cid].LatE;
        // memcpy(&map,&chain[cid].contactsLatE, sizeof(chain[cid].contactsLatE));
        memcpy(&refmap,&chref[cid].contactsLatE, sizeof(chref[cid].contactsLatE));
        break;
    case 3:
        inter_chain = chref[cid].LonS;
        // memcpy(&map,&chain[cid].contactsLonS, sizeof(chain[cid].contactsLonS));
        // printf("map 1>> \n");
        // print_contact_map(chain[cid].contactsLonS);
        // print_contact_map(map);
        memcpy(&refmap,&chref[cid].contactsLonS, sizeof(chref[cid].contactsLonS));
        // printf("refmap 1>> \n");
        // print_contact_map(chref[cid].contactsLonS);
        // print_contact_map(refmap);
        // memcpy(map, chain[cid].contactsLonS, sizeof(chain[cid].contactsLonS));
        // memcpy(refmap, chref[cid].contactsLonS, sizeof(chref[cid].contactsLonS));
        break;
    case 4:
        inter_chain = chref[cid].LatW;
        // memcpy(&map,&chain[cid].contactsLatW, sizeof(chain[cid].contactsLatW));
        memcpy(&refmap,&chref[cid].contactsLatW, sizeof(chref[cid].contactsLatW));
        break;
    default:
        std::cout << "Invalid Selection. Please try Again." << std::endl;
    } //switch
    if ( inter_chain == -1 ) {
        // exit(1);
        debug("no interface found! (%d)\n",sel);
        return;
    }

    // printf("printing maps .....................\n");
    // print_contact_map(map);
    // print_contact_map(refmap);

    // for ( int i=0; i<chref.num_atoms_ca; i++ ) {


    // class Contact {
    // public:
    //     /* int resid; */
    //     int cresid;

    //     int index; // index
    //     int cindex; // contacted index

    //     double eh; // energy scaling (well depth)
    //     double distance;


    // void print_contact_map(Contact (*map)[MAX_CONTACTS]) {
    //     int count = 0;
    //     for( int i=0; i<MOLECULE_INITIAL_SIZE; i++) {
    //         for ( int j=0; j<MAX_CONTACTS; j++ ) {
    //             if ( map[i][j].cresid != -1 ) {
    //                 printf("contact: %d %d %d %f   eh:%f\n",i,j,map[i][j].cresid,\
    //                        map[i][j].distance,map[i][j].eh);
    //                 // printf("contact: %d %d %d %f\n",i,j,map[i][j].cresid,map[i][j].distance);
    //                 count ++;
    //             }
    //         }
    //     }
    //     printf("total_contacts_printed: %d\n",count);
    // }

    debug("%d\n",chref->num_atoms_ca);

    // for (int i=0; i<chref->num_atoms_ca; i++ ) {
    for (int i=0; i<MOLECULE_INITIAL_SIZE; i++ ) {
        for ( int j=0; j<MAX_CONTACTS; j++ ) {
            if ( refmap[i][j].cresid != -1) {
                debug("contactd: %d %d %d %f   eh:%f\n",i,j,refmap[i][j].cresid,\
                       refmap[i][j].distance,refmap[i][j].eh);
                // debug("dd: %d   %d    \n%f %f %f \n%f %f %f\n",i,refmap[i][j].cresid,\
                    debug("dd: %d   %d    \n%f %f %f \n",i,refmap[i][j].cresid, \
                    chref->pos[i].x,chref->pos[i].y,chref->pos[i].z);
                //     chref->pos[refmap[i][j].cresid]);
                match_count ++;
            }
        }
        // debug("match_count: %d\n",match_count);
    // }
        // print_contact_map(refmap);
    // print_contact_map(chref->contactsLatE);
    // print_contact_map(chref->contactsLonS);
    // print_contact_map(chref->contactsLatW);
    }
    debug("end of original contacts .. total_match_count: <%d>\n",match_count);

    // for ( int j=0; j<MAX_CONTACTS; j++ ) {
    //     if ( map[i][j].cresid == -1 ) {
    //             continue;
    //         }
    //         // printf("%d %f\n",map[i][j].cresid,map[i][j].distance);
    //         // int atom = i;
    //         int first = i;
    //         int second = map[i][j].cresid;

    //         for ( int k=0; k<MAX_CONTACTS; k++ ) {
    //             if ( refmap[first][k].cresid == second ) {
    //                 match_count++;
    //                 break;
    //             }
    //             if ( refmap[second][k].cresid == first ) {
    //                 match_count++;
    //                 break;
    //             }
    //         } // for3
    //     } // for2
    // } // for1

    // PRINT HERE
    // debug("printing reference map\n");
    // print_contact_map(refmap);
    // debug("printing contact map\n");
    // print_contact_map(map);
    // debug("match_total: %d\n",match_count);


    switch(sel) {
    case 1:
        chain[cid].LonN_total2x = match_count;
        break;
    case 2:
        chain[cid].LatE_total2x = match_count;
        break;
    case 3:
        chain[cid].LonS_total2x = match_count;
        break;
    case 4:
        chain[cid].LatW_total2x = match_count;
        break;
    default:
        std::cout << "Invalid Selection. Please try Again." << std::endl;
    } //switch
}



void fprintf_contacts_for_one_monomer_face(Chain *chain, int num_chains, int sel) {

    int c = 0;

    FILE * fp_inter_contacts;

    switch(sel) {
    case 1:
        fp_inter_contacts = fopen ("contact_latlon_n.dat", "a+");
        break;
    case 2:
        fp_inter_contacts = fopen ("contact_latlon_e.dat", "a+");
        break;
    case 3:
        fp_inter_contacts = fopen ("contact_latlon_s.dat", "a+");
        break;
    case 4:
        fp_inter_contacts = fopen ("contact_latlon_w.dat", "a+");
        break;
    default:
        std::cout << "Invalid Selection. Please try Again." << std::endl;
    } //switch


    for (int i=0; i<num_chains; i++ ) {
        switch(sel) {
        case 1:
            c = chain[i].LonN_total2x;
            break;
        case 2:
            c = chain[i].LatE_total2x;
            break;
        case 3:
            c = chain[i].LonS_total2x;
            break;
        case 4:
            c = chain[i].LatW_total2x;
            break;
        default:
            std::cout << "Invalid Selection. Please try Again." << std::endl;
        } //switch


        // FPRINTF
        // debug(" %4d\n",c);
        fprintf(fp_inter_contacts," %4d",c);
    }
    debug("<%d>\n",sel);
    fprintf(fp_inter_contacts,"\n");
    fclose(fp_inter_contacts);
}


// void assign

// void build_latlon_contact_map(Chain *chain, int num_chains,char *sel) {
// void build_latlon_contact_map(Chain *chain, int num_chains, \
//                               Contact * (*map)[MOLECULE_INITIAL_SIZE][MAX_CONTACTS],int sel) {
// void build_latlon_contact_map(Chain *chain, int num_chains,Contact (*map)[MAX_CONTACTS], \
//                               Contact (*mapinter)[MAX_CONTACTS],int sel) {
void build_latlon_contact_map(Chain *chain, int num_chains, int start, int stop, \
                              Contact (*map)[MAX_CONTACTS],int sel) {

    // printf("map: %d\n",map[10]->cresid);

    // N E S W | 1 2 3 4
    // int sel = 1;
    int inter_chain;
    double dist;
    // int count = 0;

    for ( int h=start; h <stop; h++ ) {

        switch(sel) {
        case 1:
            inter_chain = chain[h].LonN;
            printf("selected LonN\n");
            // inter_chain = chain->LonN;
            break;
        case 2:
            inter_chain = chain[h].LatE;
            printf("selected LonN\n");
            // inter_chain = chain->LatE;
            break;
        case 3:
            inter_chain = chain[h].LonS;
            printf("selected LonN\n");
            // inter_chain = chain->LonS;
            break;
        case 4:
            inter_chain = chain[h].LatW;
            printf("selected LonN\n");
            // inter_chain = chain->LatW;
            break;
        default:
            std::cout << "Invalid Selection. Please try Again." << std::endl;
        } //switch
    }




    // if ( inter_chain == -1 ) {
    //     continue;
        // }
        // // printf("entered build_latlon_contact_map: %d\n",chain[h].chainid);
        // // printf("evaluating %d %d %d\n",h,inter_chain,sel);


        // for(int i=0; i < chain[h].num_atoms_ca; i++ ) {
        //     printf(" %d ",i);
        //     // need to go through whole thing .. !!!!!! ?????
        //     // or int j = i
        //     for(int j = 0; j < chain[inter_chain].num_atoms_ca; j++) {
        //         printf(" %d ",j);

        //         if ( i >= j-2 && i <= j+2) {
        //             // exclude j = i +/- 2
        //         /* printf("%d\t no --> %d j: %d %d\n",i,j,j-2,j+2); */
        //         } else {
        //         /* printf("%d\n",j); */
        //             dist = distance( chain[h].pos[i], chain[inter_chain].pos[j]);
        //             // printf("%8.4f\n",dist);

        //             if ( dist <= 8.0 ) {
        //                 /* printf("criterion-> %d %d %8.4f\n",i,j,dist); */

        //                 /* SINGLE */
        //                 /* for( count = 0; count < 40; count++ ) { */
        //                 /*     if ( mol->contact_map[i][count].cresid != -1 ) { */
        //                 /*         continue; */
        //                 /*     } else { */
        //                 /*         mol->contact_map[i][count].cresid = j; */
        //                 /*         mol->contact_map[i][count].distance = dist; */
        //                 /*         break; */
        //                 /*     } */
        //                 // }

        //                 /* DOUBLE */
        //                 for( int count=0; count < MAX_CONTACTS; count++ ) {
        //                     // if (map->[count]->cresid != -1 ) {
        //                     //     continue;
        //                     // } else {
        //                     //     map[count]->cresid = j;
        //                     //     break;
        //                     // }
        //                 }
        //                 /* and j */
        //                 // for( count = 0; count < MAX_CONTACTS; count++ ) {

        //                 //     if (mol->contacts[j][count].cresid != -1 ) {
        //                 //         continue;
        //                 //     } else {
        //                 //         mol->contacts[j][count].cresid = i;
        //                 //         mol->contacts[j][count].distance = dist;
        //                 //         break;
        //                 //     }
        //                 // }


        //                 /* break; */
        //                 /* mol->topology[i][ */
        //                 /* } else { */
        //                 /*     count */
        //                 /* } */
        //                 /* } */
        //             } /* END 8.0 */
        //         } // if i,i+1,i+2 excluded
        //     } // j-calpha inter
        // } // i-calpha
    // } // h all chains
    // /* printf("the count was: %d\n", count); */
}

// void get_8_nearestneighbors_by_centroid_from_centroids( Chain *chain, Vector arr_centroid[] ) {
//     int i;

//     // int a[7];
//     // std::cout << "Length of array = " << (sizeof(arr_centroid)/sizeof(*arr_centroid)) << std::endl;

//     // for ( i=0; i <)


// }


// ------------------PRINT HERE-----------------------------------------------
void print_contact_map(Contact (*map)[MAX_CONTACTS]) {
    int count = 0;
    for( int i=0; i<MOLECULE_INITIAL_SIZE; i++) {
        for ( int j=0; j<MAX_CONTACTS; j++ ) {
            if ( map[i][j].cresid != -1 ) {
                printf("contact: %d %d %d %f   eh:%f\n",i,j,map[i][j].cresid,\
                       map[i][j].distance,map[i][j].eh);
                // printf("contact: %d %d %d %f\n",i,j,map[i][j].cresid,map[i][j].distance);
                count ++;
            }
        }
    }
    printf("total_contacts_printed: %d\n",count);
}

// void assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh) {

void assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh,int low_resid, \
                       int high_resid,int low_cresid,int high_cresid) {
    int count = 0;
    for( int i=0; i<MOLECULE_INITIAL_SIZE; i++) {
        for ( int j=0; j<MAX_CONTACTS; j++ ) {
            if ( map[i][j].cresid != -1 ) {
                // debug("contact: %d %d %d %f  ||  eh:%f",i,j,map[i][j].cresid,\
                //        map[i][j].distance,map[i][j].eh);

                if ( i >= low_resid && i <= high_resid && \
                     map[i][j].cresid >= low_cresid &&
                    map[i][j].cresid <= high_cresid) {

                    map[i][j].eh = eh;
                    count ++;
                }
            }
        }
    }
    debug("total contact eh's changed(%d) in range %d-%d .. %d-%d\n",\
          count,low_resid,high_resid,low_cresid,high_cresid);
}
void assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh) {
    int count = 0;
    for( int i=0; i<MOLECULE_INITIAL_SIZE; i++) {
        for ( int j=0; j<MAX_CONTACTS; j++ ) {
            if ( map[i][j].cresid != -1 && map[i][j].eh == 0.0 ) {
                // debug("contact: %d %d %d %f  ||  eh:%f",i,j,map[i][j].cresid,\
                //       map[i][j].distance,map[i][j].eh);

                    map[i][j].eh = eh;
                    count ++;
            }
        }
    }
    debug("total contact eh's changed! (%d)\n",count);
}


// void fprintf_contact_map_gsop(Chain *chain,int cid1,int cid2,Contact (*map)[MAX_CONTACTS]) {
// void fprintf_contact_map_gsop(Contact (*map)[MAX_CONTACTS],float eh) {
void fprintf_contact_map_gsop(Contact (*map)[MAX_CONTACTS]) {
    // int count = 0;
    // for( int i=0; i<MOLECULE_INITIAL_SIZE; i++) {
    //     for ( int j=0; j<MAX_CONTACTS; j++ ) {
    //         if ( map[i][j].cresid != -1 ) {
    //             // printf("contact: %d %d %d %f\n",i,j,map[i][j].cresid,map[i][j].distance);
    //             printf("contact: %d %f\n",map[i][j].cresid,map[i][j].distance);
    //             count ++;
    //         }
    //     }
    // }
    // printf("total_contacts_printed: %d\n",count);


    // PRINT HERE
    // for( int i=0; i<chain[cid1].num_atoms_ca; i++ ) {
    //     printf("indices: %d",chain[cid1].indices[i]);
    //     printf(" [index %d findex %d]",chain[cid1].index,chain[cid1].findex);
    //     printf("  resid: %d",chain[cid1].resid[i]);
    //     printf("  actual_index: %d\n",chain[cid1].resid[i] - chain[cid1].resid[0]);
    // }


    // FPRINTF
    FILE * fp_gsop_top;
    // char filename[50];
    // chain[cid1].filename
    fp_gsop_top = fopen("topology.top","a+");

    int count = 0;
    // printf("  map: %d %d\n",map)
    for( int i=0; i<MOLECULE_INITIAL_SIZE; i++) {
        for ( int j=0; j<MAX_CONTACTS; j++ ) {
            if ( map[i][j].cresid != -1 ) {
                if (map[i][j].index <= map[i][j].cindex) {
                    // printf("  %3d   %3d     1    %1.5f %1.5f\n",\
                    //        map[i][j].index,map[i][j].cindex,map[i][j].distance);
                    // debug("  %3d   %3d     1    %1.5f    %1.5f",\
                    //       map[i][j].index,map[i][j].cindex,map[i][j].distance,eh);
                    debug("  %3d   %3d     1    %1.5f    %1.5f\n",\
                          map[i][j].index,map[i][j].cindex,map[i][j].distance,\
                          map[i][j].eh);

                    // fprintf(fp_gsop_top,"  %3d   %3d     1    %1.5f    %1.5f\n",\
                    //         map[i][j].index,map[i][j].cindex,map[i][j].distance,eh);
                    fprintf(fp_gsop_top," %6d  %6d    1    %1.5f    %1.5f\n",\
                            map[i][j].index,map[i][j].cindex,map[i][j].distance,\
                        map[i][j].eh);

                    count++;
                }
            }
        }
    }
    // fprintf(fp_gsop_top,"\n");
    fclose(fp_gsop_top);
    printf("total_contacts_printed: %d\n",count);
}

void fprintf_all_indices(Chain *chain, int num_chains) {

    // FILE * fp_;
    // fp_ = fopen ("file.txt", "w+");
    // fprintf(fp_, "%s %s %s %d", "We", "are", "in", 2012);
    // fclose(fp_);

    // FILE * fp_indices;
    // std::string filename = "indices_foreach_segment_%d";
    // filename
    // std::string filename << "indices_foreach_segment_"
    //                      << num_chains << ".dat";
    // printf("fn: %s",filename);


    //  // Taking the string value :
    // int n_chain = num_chains - 1;
    // std::string filename;
    // filename = str( boost::format("indices_foreach_segment_0_%d.dat") % n_chain);
    // // std::cout << filename << "\n" << std::endl;
    // // printf("filename: %s\n",filename.c_str());


    // // fp_indices = fopen(filename.c_str(),"w+");
    // // for ( int i=0; i < num_chains; i++ ) {
    // //     fprintf(fp_indices,"%6d  %6d  %6d\n",\
    // //             chain[i].chainid,chain[i].index,chain[i].findex);
    // // }
    // // fclose(fp_indices);


    // // #include <fstream>
    // std::ofstream myfile;
    // myfile.open("indices_%d",n_chain;
    // // std::ofstream fout;



    // for ( int i=0; i < num_chains; i++ ) {
    //     // fprintf(fp_indices,"%6d  %6d  %6d\n",\
    //     //         chain[i].chainid,chain[i].index,chain[i].findex);

    //     myfile << std::setw(6) << chain[i].chainid
    //            << std::setw(6) << chain[i].index
    //            << std::setw(6) << chain[i].findex
    //            << std::endl;
    // }

    // myfile.close();


    std::ostringstream fname;
    // fname << std::setw(3) << std::setfill('0') << 5 << ".dat";
    fname << "indices_foreach_segment_0_" << num_chains - 1 << ".dat";
    std::ofstream fout(fname.str().c_str(), std::ios::out);

    for ( int i=0; i < num_chains; i++ ) {
        // fprintf(fp_indices,"%6d  %6d  %6d\n",\
        //         chain[i].chainid,chain[i].index,chain[i].findex);

        printf("%6d  %6d  %6d\n",                                   \
               chain[i].chainid,chain[i].index,chain[i].findex);

        fout << std::setw(6) << chain[i].chainid
             << std::setw(6) << chain[i].index
             << std::setw(6) << chain[i].findex
             << std::endl;
    }
    fout.close();
}



// void print_contact_map_gsop(Contact (*map)[MAX_CONTACTS]) {
//     int count = 0;
//     for( int i=0; i<MOLECULE_INITIAL_SIZE; i++) {
//         for ( int j=0; j<MAX_CONTACTS; j++ ) {
//             if ( map[i][j].cresid != -1 ) && (map[i][j].cresid ) {
//                 printf("contact: %d %d %d %f\n",i,j,map[i][j].cresid,map[i][j].distance);
//                 count ++;
//             }
//         }
//     }
//     printf("total_contacts_printed: %d\n",count);
// }

// void load_dcd_coords_to_chain(dcdhandle *dcd,Chain *chain,int num_chains) {

//     for (int i=0; i<num_chains; i++) {

//         for ( int j=0; j<chain[i].num_atoms_ca; j++ ) {
//             // printf("%d ",chain[i].indices[j]);
//             chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
//             chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
//             chain[i].pos[j].z = dcd->z[chain[i].indices[j]];
//             // chain[i].pos[j].y = double (dcd->y[chain[i].indices[j]]);
//             // chain[i].pos[j].z = double (dcd->z[chain[i].indices[j]]);

//             // chain[i].pos[chain[i].indices[j]].y = double (dcd->y);
//             // chain[i].pos[chain[i].indices[j]].z = double (dcd->z);
//             // chain->pos[i].x = x;
//             // chain->pos[i].y = y;
//             // chain->pos[i].z = z;
//         }
//     }
// }

// void load_chain_coords_to_timestep(dcdhandle *dcd,Chain *chain,int num_chains){
// void load_chain_coords_to_timestep(dcdhandle *dcd,Chain *chain,int num_chains, \
//                                    void *v, const molfile_timestep_t *ts,int natoms){
// void load_chain_coords_to_timestep(Chain *chain,int num_chains,const molfile_timestep_t *ts){
//     debug("the number of chains to load: %d\n",num_chains);
//     debug("maximum_index: %d\n",chain[num_chains-1].findex);
//     debug("chain_count: %d\n",num_chains);

//     // for(int i=0; i<num_chains; i++){
//     //     for(int j=0; j<chain[i].num_atoms_ca; j++){
//     //         printf("%f %f %f\n",chain[i].pos[j].x,\
//     //                chain[i].pos[j].y,\
//     //                chain[i].pos[j].z);
//     //     }
//     //     printf("next_chain!\n");
//     // }
//     // printf("done with chain.");
//     // exit(0);
//     // return;

//     // Works! - not necessary
//     // dcdhandle *dcd1;
//     // dcd1 = (dcdhandle *)malloc(sizeof(dcdhandle));
//     // memset(dcd1, 0, sizeof(dcdhandle));
//     // dcd1->x = (float *)malloc(natoms * sizeof(float));
//     // dcd1->y = (float *)malloc(natoms * sizeof(float));
//     // dcd1->z = (float *)malloc(natoms * sizeof(float));

//     // int x,c,sindex,findex,total,natoms;
//     int natoms,count,count1;
//     natoms = count = count1 = 0;
//     // x = c = sindex = findex = total = natoms = 0;
//     // c will be current chain. 0,1

//     for(int h=0; h<num_chains; h++){
//         // debug("chain: %d\n",h);
//         // debug("atoms in chain(ca) %d\n",chain[h].num_atoms_ca);
//         natoms += chain[h].num_atoms_ca;
//         // sindex = chain[h].index;
//         // findex = chain[h].findex;

//         for(int i=0; i<chain[h].num_atoms_ca; i++){
//             // OLD
//             // printf("c %d %d %d\n",count,count+1,count+2);
//             // printf("i,sindex,findex: %d %d %d\n",i,sindex,findex);
//             // printf("%f ",chain[c].pos[i].x);
//             // printf("%f ",chain[c].pos[i].y);
//             // printf("%f \n",chain[c].pos[i].z);
//             // ts->coords[x] = (float)chain[c].pos[sindex].x;
//             // ts->coords[x+1] = (float)chain[c].pos[sindex].y;
//             // ts->coords[x+2] = (float)chain[c].pos[sindex].z;

//             // printf("%f")


//             // LOADING
//             // debug("xyz: %f %f %f\n",chain[h].pos[i].x,  \
//             //        chain[h].pos[i].y,\
//             //        chain[h].pos[i].z);

//             ts->coords[count] = (float)chain[h].pos[i].x;
//             ts->coords[count+1] = (float)chain[h].pos[i].y;
//             ts->coords[count+2] = (float)chain[h].pos[i].z;
//             count += 3;
//             count1 ++;
//             // sindex++;
//         }
//     }
//     debug("number of atoms: %d\n",natoms);
//     debug("total coordinates written: %d  %d\n",count1,count);

//     // CHECK TS
//     // for(int i=0; i<natoms; i++){
//     //     printf("i: %d   %f %f %f\n",i*3,ts->coords[i*3],ts->coords[i*3+1],ts->coords[i*3+2]);
//     // }

//     return;
// }

void print_chain_coords(Chain *chain,int num_chains) {
    for (int i=0; i<num_chains; i++) {
        // printf("chain:%d\n",i);
        for ( int j=0; j<chain[i].num_atoms_ca; j++ ) {

            // PRINT HERE
            printf("index: %d >",chain[i].indices[j]);
            printf("xyz: %f %f %f\n",chain[i].pos[j].x,\
                   chain[i].pos[j].y,chain[i].pos[j].z);
            break;


            // chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
            // chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
            // chain[i].pos[j].z = dcd->z[chain[i].indices[j]];
            // chain[i].pos[j].y = double (dcd->y[chain[i].indices[j]]);
            // chain[i].pos[j].z = double (dcd->z[chain[i].indices[j]]);

            // chain[i].pos[chain[i].indices[j]].y = double (dcd->y);
            // chain[i].pos[chain[i].indices[j]].z = double (dcd->z);
            // chain->pos[i].x = x;
            // chain->pos[i].y = y;
            // chain->pos[i].z = z;
        }
    }
}
