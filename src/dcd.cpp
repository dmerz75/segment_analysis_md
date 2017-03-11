// dcd.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h> // strtod?, stod


// namespaces
// using namespace

// extern
extern int endoffile;

/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
// #include "readfile.h"
#include "chain.h"
#include "md.h"
/* #include "dcdio.h" */


/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
int advance_dcd(int numframes,int frame,dcdhandle *v,int natoms,molfile_timestep_t *timestep) {
    for (int i=0; i<numframes; i++) {
        int rc = read_next_timestep(v, natoms, timestep);
        if (rc) {
            fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
            break;
            // exit(1);
        }

        if ( i == frame) {
            break;
        }
    }
    return frame + 1;
    // return frame;
}
void load_chain_to_timestep(Chain *chain,int num_chains,const molfile_timestep_t *ts){
    debug("the number of chains to load: %d\n",num_chains);
    debug("maximum_index: %d\n",chain[num_chains-1].findex);
    debug("chain_count: %d\n",num_chains);

    // Works! - not necessary
    // dcdhandle *dcd1;
    // dcd1 = (dcdhandle *)malloc(sizeof(dcdhandle));
    // memset(dcd1, 0, sizeof(dcdhandle));
    // dcd1->x = (float *)malloc(natoms * sizeof(float));
    // dcd1->y = (float *)malloc(natoms * sizeof(float));
    // dcd1->z = (float *)malloc(natoms * sizeof(float));

    // int x,c,sindex,findex,total,natoms;
    int natoms,count,count1;
    natoms = count = count1 = 0;

    for(int h=0; h<num_chains; h++){

        natoms += chain[h].num_atoms_ca;

        for(int i=0; i<chain[h].num_atoms_ca; i++){
            ts->coords[count] = (float)chain[h].pos[i].x;
            ts->coords[count+1] = (float)chain[h].pos[i].y;
            ts->coords[count+2] = (float)chain[h].pos[i].z;
            count += 3;
            count1 ++;
        }
    }
    debug("number of atoms: %d\n",natoms);
    debug("total coordinates written: %d  %d\n",count1,count);

    // CHECK TS
    // for(int i=0; i<natoms; i++){
    //     printf("i: %d   %f %f %f\n",i*3,ts->coords[i*3],ts->coords[i*3+1],ts->coords[i*3+2]);
    // }

    return;
}

void load_dcd_to_chain(dcdhandle *dcd,Chain *chain,int num_chains) {

    for (int i=0; i<num_chains; i++) {

        for ( int j=0; j<chain[i].num_atoms_ca; j++ ) {
            // printf("%d ",chain[i].indices[j]);
            chain[i].pos[j].x = dcd->x[chain[i].indices[j]];
            chain[i].pos[j].y = dcd->y[chain[i].indices[j]];
            chain[i].pos[j].z = dcd->z[chain[i].indices[j]];
        }
    }
}
