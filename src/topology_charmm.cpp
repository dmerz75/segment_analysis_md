// topology_charmm.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
// #include <stdio.h>
// #include <stdlib.h> // strtod?, stod
// #include <assert.h>
// #include <string.h>
// #include <string>
// #include <iostream>
// #include <string>


// #include <fstream>
// #include <stdlib.h>
// #include <stdio.h>
// #include <algorithm> // count
// #include <iterator> // istream_iterator
// #include <ctime>

// namespaces
// using namespace

// extern
// extern int endoffile;

/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
// #include "readfile.h"
#include "chain.h"
// #include "md.h"

// Vector pos;

// #define MOLECULE_INITIAL_SIZE 700
// Chain chain;


/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
void top_charmm_write_general(Chain *chain,int max_num_chains,FILE *fp) {
    printf("Hello from Topology Charmm!\n");

    fprintf(fp,"MASS %d");
    // MASS    20 C     12.01100 C ! carbonyl C, peptide backbone

    // resname
    for ( int c=0; c < max_num_chains; c++ ) {
        for ( int i=0; i < chain[c].num_atoms; i++ ) {
            printf("%d ",i);
            printf("%s ",chain[c].resname[i]);
            printf("%d \n",chain[c].resid[i]);
            fprintf(fp,"%d %s %d\n",i,chain[c].resname[i],chain[c].resid[i]);
        }
        printf("next chain.\n");
    }
}
