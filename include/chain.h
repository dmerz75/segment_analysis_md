#ifndef _CHAIN_H_
#define _CHAIN_H_



/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>

// my headers
#include <md.h>
#include <debug.h>


#include "dcdio.h"


// headers 3rd
/* #include <armadillo> */
/* #include <gsl/gsl_math.h> */
/* #include <gsl/gsl_eigen.h> */

/* #include "boost/multi_array.hpp" */
#include <cassert>

/* SYMLINK: Eigen /usr/local/include/. */
/* [dale~/sop_dev/contacts/segment]$lt /usr/local/include/eigen                        */
/* lrwxrwxrwx 1 root root 25 09.25.2015 21:04 /usr/local/include/eigen -> /usr/include/eigen3/Eigen/ */
/* #include<Eigen/Eigenvalues> */
/* #include <Eigen/Dense> */

#ifdef INERTIA
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
using namespace Eigen;
#endif // INERTIA

/* using namespace Eigen; */



 // namespaces
/* using namespace arma; */

/* ---------------------------------------------------------
   definitions
   --------------------------------------------------------- */
/* 2000, 100;  700,20; 500,50 */
#if defined ONEMOL || defined ALLATOM
#define MOLECULE_INITIAL_SIZE 6900
#define MOLECULE_INITIAL_SIZE_MT 6900
#else
#define MOLECULE_INITIAL_SIZE 700
#define MOLECULE_INITIAL_SIZE_MT 700
#endif // ONEMOL
#define MAX_CONTACTS 24 /* per atom */

class Contact {
public:
    /* int resid; */
    int cresid;

    int total;

    int index; // index
    int cindex; // contacted index

    double eh; // energy scaling (well depth)
    double distance;


    // Constructor
    Contact (); // Constructor declared.


    // Destructors:
    ~Contact() {
        /* delete [] pos; */
    };

};
inline Contact::Contact() {
    cresid = -1;
    index = -1;
    cindex = -1;
    total = -1;
    eh = 0.0; // keep at 0.0, used in assign_contact_eh as default
    distance = 0.0;
}



typedef struct Chain Molecule;


class Chain {
public:
    char filename[46];

    int chainid;
    int num_atoms;
    int num_atoms_ca;

    int index;
    int findex;
    int pf_pos_low; // for selecting index near midpoint of centroids, curvature of PF
    int pf_pos_high;
    int file_line_begin;
    int file_line_end;

    char atomtype[MOLECULE_INITIAL_SIZE][4];
    char resname[MOLECULE_INITIAL_SIZE][4];
    int resid[MOLECULE_INITIAL_SIZE];
    int indices[MOLECULE_INITIAL_SIZE]; // not used yet.

    Vector pos[MOLECULE_INITIAL_SIZE];
    Contact contacts[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];

    int total_num_contacts_2x;
    int tally_num_contacts_2x[MOLECULE_INITIAL_SIZE];
    int num_contacts_2x_persist_ref[MOLECULE_INITIAL_SIZE];

    Vector centroid;
    int nearest_targets[32];
    double nearest_targets_dist[32];
    int nearestneighbors_chainid[8];
    double nearestneighbors_dist[8];


    int LonN; // contacted monomer
    int LatE;
    int LonS;
    int LatW;
    int LonN_total2x; // total
    int LatE_total2x;
    int LonS_total2x;
    int LatW_total2x;


    // general inter
    Contact intercontacts[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];

    // latlon contacts
    Contact contactsLonN[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    Contact contactsLatE[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    Contact contactsLonS[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    Contact contactsLatW[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];



    // Constructors:
    Chain(int a=-1, int b = 0, double p = 0.0, double q = -1.0) {
        chainid = a;
        num_atoms = a;
        num_atoms_ca = a;
        index = a;
        findex = a;
        pf_pos_low = a;
        pf_pos_high = a;
        file_line_begin = a;
        file_line_end = a;
        total_num_contacts_2x = a;

        centroid.x = p;
        centroid.y = p;
        centroid.z = p;

        LonN = a;
        LonS = a;
        LatW = a;
        LatE = a;

        // initialize to 0.0 p
        for (int i=0; i < MOLECULE_INITIAL_SIZE; i++ ) {
            pos[i].x = p;
            pos[i].y = p;
            pos[i].z = p;

            // initialize to 0 b
            tally_num_contacts_2x[i] = b;
            num_contacts_2x_persist_ref[i] = b;


            for ( int j=0; j < MAX_CONTACTS; j++ ) {
                contacts[i][j].cresid = a;
                contacts[i][j].distance = p;
            }
            for ( int j=0; j < 8; j++ ) {
                nearestneighbors_chainid[j] = -1;
                nearestneighbors_dist[j] = p;
            }
            for ( int j=0; j < 32; j++ ) {
                nearest_targets[j] = a;
            }
        }
    }

    // Destructors:
    ~Chain() {
        /* delete [] pos; */
        /* std::cout << "Object is being deleted" << std::endl; */
    };

    // Declarations:
    void print_prop();
    void assign_indices();
    void fprintf_indices();
    void GetMomentofInertiaTensor();
    void ComputeCentroid();
    void print_centroid();

    int countcontacts(Contact (*map)[MAX_CONTACTS]);

};
inline void Chain::ComputeCentroid() {
    int i;
    double xtot = 0.0, ytot = 0.0, ztot = 0.0;
    /* printf("computing centroid ca: %d\n",num_atoms_ca); */
    for ( int i=0; i < num_atoms_ca; i++ ) {
        /* printf("%d\n",i); */
        xtot += pos[i].x;
        ytot += pos[i].y;
        ztot += pos[i].z;
    }
    centroid.x = xtot / num_atoms_ca;
    centroid.y = ytot / num_atoms_ca;
    centroid.z = ztot / num_atoms_ca;
    // PRINT HERE
    /* debug("c(%d): %f %f %f\n",chainid,centroid.x,centroid.y,centroid.z); */
}

// amazing!!
/* By declared that method with inline keyword, compiler will either inline */
/*     the whole method or, if it decides not to, it will generate anonymous */
/*     method (same method but with some unique name for a given object file), */
/*     so there will be no conflicts in object files. For example: */
/* class foo { */
/* public: */
/*   void bar (); */
/* }; */
/* inline void foo::bar () */
/* { */
/* } */
inline void Chain::print_prop() {
    /* printf("id: %d ",chainid); */
    /* printf(" filename: %s\n",filename); */
    /* printf("CA: %d ",num_atoms_ca); */
    /* printf(" index(%d) ",index); */
    /* printf(" final_index(%d)\n",findex); */
    /* printf("begin_end: [ %d %d ]\n\n",file_line_begin,file_line_end); */

    /* debug("id: %d   CA: %d\n",chainid,num_atoms_ca); */
    /* debug(" filename: %s\n",filename); */
    /* debug("CA: %d \n",num_atoms_ca); */
    /* debug(" index(%d) \n",index); */
    /* debug(" final_index(%d)\n",findex); */
    /* debug("begin_end: [ %d %d ]\n",file_line_begin,file_line_end); */

    printf("centroid: %f %f %f\n",centroid.x,centroid.y,centroid.z);
}
inline void Chain::print_centroid() {
    printf("centroid: %f %f %f\n",centroid.x,centroid.y,centroid.z);
}
inline void Chain::assign_indices() {
    /* printf("id: %d ",chainid); */
    /* printf(" filename: %s\n",filename); */
    /* printf("CA: %d ",num_atoms_ca); */
    /* printf(" index(%d) ",index); */
    /* printf(" final_index(%d)\n",findex); */
    /* printf("begin_end: [ %d %d ]\n\n",file_line_begin,file_line_end); */

    debug("id: %d   atoms:%d   CA:%d\n",chainid,num_atoms,num_atoms_ca);
    /* debug(" filename: %s\n",filename); */
    /* debug("CA: %d \n",num_atoms_ca); */
    /* debug(" index(%d) \n",index); */
    /* debug(" final_index(%d)\n",findex); */
    debug("begin_end: [ %d %d ]\n",file_line_begin,file_line_end);

    for( int i=0; i < num_atoms_ca; i++ ) {
        indices[i] = index + i;
    }
}
inline void Chain::fprintf_indices() {
    /* printf(" index(%d) ",index); */
    /* printf(" final_index(%d)\n",findex); */
    printf(" %d__[%d|%d]\n",chainid,index,findex);

    FILE * fp_indices;
    fp_indices = fopen("indices_foreach_segments.dat","a+");
    for ( int i=0; i < num_atoms_ca; i++ ) {
        fprintf(fp_indices,"%6d  %6d  %6d\n",chainid,index,findex);
    }
    fclose(fp_indices);

    /* for( int i=0; i < num_atoms_ca; i++ ) { */
    /*     /\* indices[i] = index + i; *\/ */
    /*     printf(" %d",indices[i]); */
    /*     printf(" ") */
    /* } */
    /* printf("\n"); */
}

inline int Chain::countcontacts(Contact (*map)[MAX_CONTACTS]) {

    /* printf("declare countcontacts!\n"); */
    int count = 0;

    // using 0,3 indices. 1,4, 2,5, ... 430,433
    for(int i=0; i<MOLECULE_INITIAL_SIZE; i++){
        // for(int j=0; j<MAX_CONTACTS; j++){
        for(int k=0; k<MAX_CONTACTS; k++){

            if(map[i][k].cresid == -1){
                break;
            } else {
                count ++;
                // printf("%d-%d (%d)  %d/50\n",i,map[i][k].cresid,count,k);
            }
            // printf(" %d",map[i][k].index);
            // map[num1][k].cindex = chain[num1].indices[j];
            // double eh; // energy scaling (well depth)
            // map[num1][k].distance = dist;

            // if(map[i][k].cresid != -1){
            //     count ++;
            // } // if
            // } // for 3
            // } // if
        } // for 2
    } // for 1

    if(count != 0){
        debug("declare total contacts in this map are: %d\n",count);
    } else {
        debug("empty map.\n");
    }
    // // map[chain_num].total = count;

    return count;
}

#ifdef INERTIA
inline void Chain::GetMomentofInertiaTensor() {

    /* printf("chain: %d;  atoms: %d;\n",chainid,num_atoms_ca); */
    double x2=0.0,y2=0.0,z2=0.0;
    double xy=0.0,xz=0.0,yz=0.0;

    double Ixx=0.0,Iyy=0.0,Izz=0.0;
    double Ixy=0.0,Ixz=0.0,Iyz=0.0;


    /* http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node64.html moment of Inertia Tensor */
    /* http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node67.html get eigenvectors */
    for(int a=0; a<num_atoms_ca; a++) {
    /* for(int a=0; a<5; a++) { */
        /* printf("%d",a); */
        /* printf("%f %f %f\n",pos[a].x,pos[a].y,pos[a].z); */
        x2 += pos[a].x * pos[a].x;
        y2 += pos[a].y * pos[a].y;
        z2 += pos[a].z * pos[a].z;

        xy += pos[a].x * pos[a].y;
        xz += pos[a].x * pos[a].z;
        yz += pos[a].y * pos[a].z;

        Ixx += y2 + z2;
        Iyy += x2 + z2;
        Izz += x2 + y2;

        Ixy += pos[a].x * pos[a].y;
        Ixz += pos[a].x * pos[a].z;
        Iyz += pos[a].y * pos[a].z;
    }

    Ixy = Ixy * -1.0;
    Ixz = Ixz * -1.0;
    Iyz = Iyz * -1.0;


    printf("xx_yy_zz: %f %f %f\n",Ixx,Iyy,Izz);
    printf("xy_xz_yz: %f %f %f\n",Ixy,Ixz,Iyz);


    // begin eigen library

    Matrix3d A(3,3);
    A(0,0) = Ixx;
    A(0,1) = Ixy;
    A(0,2) = Ixz;

    A(1,0) = Ixy;
    A(1,1) = Iyy;
    A(1,2) = Iyz;

    A(2,0) = Ixz;
    A(2,1) = Iyz;
    A(2,2) = Izz;
    /* std::cout << "Here is the matrix A:\n" << A << std::endl; */


    // build matrix
    double InertiaTensor[3][3] = {
        {Ixx, Ixy, Ixz},
        {Ixy, Iyy, Iyz},
        {Ixz, Iyz, Izz},
    };


    // Ixy xz yz are altered beyond this point.
    double lowest_value_diag = Ixx;
    if ( Iyy < Ixx ) {
        if ( Izz < Iyy ) {
            lowest_value_diag = Izz;
        } else {
            lowest_value_diag = Iyy;
        }
    } else if ( Izz < Ixx ) {
        lowest_value_diag = Izz;
    }

    double lowest_value_offdiag;
    if ( Ixy < 0.0 ) {
        Ixy = -1.0 * Ixy;
    }
    if ( Ixz < 0.0 ) {
        Ixz = -1.0 * Ixz;
    }
    if ( Iyz < 0.0 ) {
        Iyz = -1.0 * Iyz;
    }

    lowest_value_offdiag = Ixy;
    if ( Ixz < Ixy ) {
        if ( Iyz < Ixz ) {
            lowest_value_offdiag = Iyz;
        } else {
            lowest_value_offdiag = Ixz;
        }
    } else if ( Iyz < Ixy ) {
        lowest_value_offdiag = Iyz;
    }

    double lowest_value;
    if ( lowest_value_offdiag < lowest_value_diag ) {
        lowest_value = lowest_value_offdiag;
    } else {
        lowest_value = lowest_value_offdiag;
    }

    printf("lv: %f\n",lowest_value);
    lowest_value = lowest_value * 0.95;
    /* printf("lv: %f\n",lowest_value); */
    // Matrix entries * lowest_value;



    /* // ORIGINAL EXAMPLE */
    /* // Create a 3D array that is 3 x 4 x 2 */
    /* typedef boost::multi_array<double, 3> array_type; */
    /* typedef array_type::index index; */
    /* array_type A(boost::extents[3][4][2]); */

    /* // Assign values to the elements */
    /* int values = 0; */
    /* for(index i = 0; i != 3; ++i) */
    /*     for(index j = 0; j != 4; ++j) */
    /*         for(index k = 0; k != 2; ++k) */
    /*             A[i][j][k] = values++; */

    /* // Verify values */
    /* int verify = 0; */
    /* for(index i = 0; i != 3; ++i) */
    /*     for(index j = 0; j != 4; ++j) */
    /*         for(index k = 0; k != 2; ++k) */
    /*             assert(A[i][j][k] == verify++); */
    /* /\* return 0; *\/ */



    /* ---------------------------------------------------------
       boost! no eigenvalue/eigenvector
       --------------------------------------------------------- */
    /* // Create a 3D array that is 3 x 4 x 2 */
    /* typedef boost::multi_array<double, 2> array_type; */
    /* typedef array_type::index index; */
    /* array_type A(boost::extents[3][3]); */
    /* array_type B(boost::extents[3][3]); */


    /* // x row */
    /* A[0][0] = Ixx; */
    /* A[0][1] = Ixy; */
    /* A[0][2] = Ixz; */

    /* // y row */
    /* A[1][0] = Ixy; */
    /* A[1][1] = Iyy; */
    /* A[1][2] = Iyz; */

    /* // y row */
    /* A[2][0] = Ixz; */
    /* A[2][1] = Iyz; */
    /* A[2][2] = Izz; */



    /* /\* // Verify values *\/ */
    /* /\* int verify = 0; *\/ */
    /* /\* for(index i = 0; i != 3; i++) { *\/ */
    /* /\*     for(index j = 0; j != 3; j++) { *\/ */
    /* /\*         /\\* assert(A[i][j] == verify++); *\\/ *\/ */
    /* /\*         // reassign by lowest_value scaling *\/ */
    /* /\*         printf("A: %f\n",A[i][j]); *\/ */
    /* /\*     } *\/ */
    /* /\* } *\/ */

    /* // Verify values A */
    /* for(index i = 0; i != 3; i++) { */
    /*     for(index j = 0; j != 3; j++) { */
    /*         // reassign by lowest_value scaling */
    /*         A[i][j] = A[i][j] / lowest_value; */
    /*         printf("A/d: %f\n",A[i][j]); */
    /*     } */
    /* } */
    /* ---------------------------------------------------------
       boost! no eigenvalue/eigenvector - END
       --------------------------------------------------------- */





    /* gsl_matrix_view m = gsl_matrix_view_array (InertiaTensor, 3, 3); */
    /* gsl_vector *eval = gsl_vector_alloc (3); */
    /* gsl_matrix *evec = gsl_matrix_alloc (3, 3); */

    /* gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3); */
    /* gsl_eigen_symmv (&m.matrix, eval, evec, w); */
    /* gsl_eigen_symmv_free (w); */
    /* gsl_eigen_symmv_sort (eval,evec,GSL_EIGEN_SORT_ABS_ASC); */
    /* { */
    /*     int i; */

    /*     for (i = 0; i < 4; i++) */
    /*     { */
    /*         double eval_i */
    /*             = gsl_vector_get (eval, i); */
    /*         gsl_vector_view evec_i */
    /*             = gsl_matrix_column (evec, i); */

    /*         printf ("eigenvalue = %g\n", eval_i); */
    /*         printf ("eigenvector = \n"); */
    /*         gsl_vector_fprintf (stdout, */
    /*                             &evec_i.vector, "%g"); */
    /*     } */
    /* } */

    /* gsl_vector_free (eval); */
    /* gsl_matrix_free (evec); */

    /* return 0; */



    // eigen continued
    A = A / lowest_value;
    std::cout << "Here is the matrix A:\n" << A << std::endl;

    VectorXcd eivals = A.eigenvalues();
    std::cout << "The eigenvalues of the 3x3 matrix of ones are:" << std::endl << eivals << std::endl;

    /* VectorXcd eivecs = A.eigenvectors(); */
    /* std::cout << "The eigenvalues of the 3x3 matrix of ones are:" << std::endl << eivecs << std::endl; */

    EigenSolver<MatrixXd> es(A);
    std::cout << "The first eigenvector of the 3x3 matrix of ones is:"
              << std::endl << es.eigenvectors().col(0) << std::endl;

    std::cout << "The second eigenvector of the 3x3 matrix of ones is:"
              << std::endl << es.eigenvectors().col(1) << std::endl;

    std::cout << "The third eigenvector of the 3x3 matrix of ones is:"
              << std::endl << es.eigenvectors().col(2) << std::endl;

    /* Matrix3d m = Matrix3d::(3,3); */
    /* cout << "m =" << endl << m << endl; */

}
#endif // INERTIA




/* class Map { */
/*     /\* int num_atoms,max_contacts_per_atom; *\/ */
/* public: */

/*     /\* std::vector< std::vector< std::vector<int>>> map; *\/ */
/*     std::vector< std::vector<int> > map; */
/*     /\* std::vector *\/ */
/*     std::vector<Contact> contacts; */

/*     // Constructor */
/*     /\* Map(); *\/ */

/*     /\* declarations: *\/ */
/*     /\* void extend_map(int); *\/ */
/*     void build_map_default(Chain *chain, int cid); // use chain_ref */

/* }; */
/* /\* inline Map::extend_map(int ca){ *\/ */
/* /\*     for (int i=0; i<ca; i++) { *\/ */
/* /\*         map.push_back(i); // ADD the ROW. *\/ */
/* /\*     } *\/ */
/* /\* } *\/ */
/* inline void Map::build_map_default(Chain *chain,int cid){ */
/*     // default: 8A */




/*     /\* int count; *\/ */
/*     double dist; */
/*     Contact x; */
/*         /\* ,y; *\/ */

/*     for (int i=0; i<chain[cid].num_atoms_ca; i++) { */
/*         /\* printf("%d\n",i); *\/ */

/*         for(int j=i; j<chain[cid].num_atoms_ca; j++) { */

/*             if ( i >= j-2 && i <= j+2) { */
/*                 // exclude j = i +/- 2 */
/*                 /\* printf("%d\t no --> %d j: %d %d\n",i,j,j-2,j+2); *\/ */
/*             } else { */

/*                 /\* printf("%d\n",j); *\/ */
/*                 dist = distance(chain[cid].pos[i],chain[cid].pos[j]); */
/*                 /\* printf("%8.4f\n",dist); *\/ */

/*                 if ( dist <= 8.0 ) { */

/*                     x.cresid = j; */
/*                     x.distance = dist; */
/*                     x.index = i + chain[cid].index; */
/*                     x.cindex = j + chain[cid].index; */

/*                     contacts.insert(contacts.end(),x); */
/*                 } /\* END 8.0 *\/ */
/*             } // else/if */

/*         } // j */


/*         /\* map.push_back(contacts); *\/ */
/*         /\* contacts.clear(); *\/ */
/*     } // i */

/*     /\* std::vector<int>::iterator it; *\/ */
/*     /\* for (it=contacts.begin(); it<contacts.end(); it++) { *\/ */
/*     /\*     // std::cout << *it << ' ' << std::endl; *\/ */
/*     /\*     printf("%d ",*it); *\/ */
/*     /\* } *\/ */

/*     /\* std::vector<int>::iterator iter; // declare a general iterator *\/ */
/*     for ( int k = contacts.begin(); k < contacts.end(); k++){ */
/*         printf("%d\n",contacts[k].index); */
/*     } */

/* } */


/* inline Map::Map() { */
/*     // constructor */

/*     for (int i=0; i<num_atoms; i++ ) { */
/*         /\* printf("%d\n",i); *\/ */
/*         for (int j=0; j<max_contacts_per_atom; j++){ */
/*             /\* if ( mol->contacts[i][j].cresid != 0 || mol->contacts[i][j].distance != 0.0) { *\/ */
/*             /\* printf("%d %f\n",mol->contacts[i][j].cresid,mol->contacts[i][j].distance); *\/ */
/*             contacts[i][j].cresid = -1; */
/*             contacts[i][j].distance = 0.0; */
/*         } */
/*     } */
/* } */
/* inline void Map::get_contacts() { */

/* } */
/* inline void Map::set_values (int natoms, int max_contacts) { */
/*     num_atoms = natoms; */
/*     max_contacts_per_atom = max_contacts; */

/* } */




/* ---------------------------------------------------------
   Molecule preparations:
   --------------------------------------------------------- */
void build_contact_map(Molecule *mol);
void build_latlon_contact_map(Chain *chain, int num_chains, int start, int stop,\
                              Contact (*map)[MAX_CONTACTS],int sel);
void build_contact_map_inter_monomer(Chain *chain,int cid, Contact (*map)[MAX_CONTACTS], int sel);
void build_contact_map_inter(Chain *chain,int cid1,int cid2,Contact (*map)[MAX_CONTACTS]);


void compare_contacts(Chain *chref,Chain *chain,int cid, int sel);
void evaluate_original_contacts_now(Chain *chref,Chain *chain,int cid, int sel);

void print_contact_map(Contact (*map)[MAX_CONTACTS]);
void fprintf_count_of_contacts_2x_persist_from_reference( Molecule *m1, Molecule *m2);
void fprintf_bond_vector_angle_with_tension_vector( Molecule *m1, Molecule *m2 );
void fprintf_contacts_for_one_monomer_face(Chain *chain, int num_chains, int sel);
void fprintf_contact_map_gsop(Contact (*map)[MAX_CONTACTS]);
void fprintf_all_indices(Chain *chref,int num_chains);


void tubulin_monomers_angles( Molecule *m1, Molecule *m2, int r1, int r2, FILE *fp);
void get_tubulin_centroid( Molecule *m1, Vector *com);
void get_centroid_from_array( Vector arr_centroid[],int i, int j, Vector *centroid);
Vector get_centroid_from_centroids(Chain *chain,int start_chain,int stop_chain);
void get_centroid_movement(Chain *chain_ref,Chain *chain_later,int num_chains);
void get_centroid_xyzposition_movement(Chain *chain_later,int chains_to_use,FILE *fp_xyzpos);
/* void get_centroid_xyzposition_movement(chain_later,chains_to_use,fp_xyzpos); */


/* void get_latangle_from35centroids(Chain *chain_ref,Chain *chain_later,int m,int sel); */
void get_latangle_from35centroids(Chain *chain_ref,Chain *chain_later,int m,FILE *fp,int sel);
std::vector<double> get_3angles_dimerbyalpha(Chain *chain_ref,Chain *chain_later,int cid);

/* void get_latangle_from35centroids(Chain *chain_ref,Chain *chain_later,int m,int sel,FILE *fp_ang); */


void get_midpoint_reference_indices(Chain *chain,int i,int *index_ref_i,int *index_ref_c,int sel);
/* get_midpoint_reference_indices(chain_ref,i,&index_ref_i,&index_ref_c,1); */


void compute_nearest_neighbors(Chain *chain, int num_chains);
void determine_latlon_neighbors(Chain *chains, int num_chains, Vector axis);
void get_monomer_monomer_info( Molecule *m1, Molecule *m2,\
                               Vector ref_norm, Vector m1centroid,\
                               Vector m2centroid );



/* void load_dcd_coords_to_chain(dcdhandle *dcd,Chain *chain,int num_chains); */
/* void load_chain_coords_to_timestep(Chain *chain,int num_chains,\ */
/*                                    const molfile_timestep_t *ts); */

void print_chain_coords(Chain *chain,int num_chains);



// overloaded.
void assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh,int low_resid,\
                       int high_resid,int low_cresid,int high_cresid);
void assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh);



/* void dcd_to_chain(ddc,chain_0); */
/* void dcd_to_chain(dcdhandle *v,Chain *chain); */

/* void load_chain_coords_to_timestep(dcdhandle *dcd,Chain *chain,int num_chains,\ */
/*                                    void *v, const molfile_timestep_t *ts,int natoms); */
/* load_chain_coords_to_timestep(dcdw,chain_later,chains_to_use,vw,&timestep_w,natoms_w); */

/* static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) { */
/* static int write_timestep(void *v, const molfile_timestep_t *ts) { */





#endif
