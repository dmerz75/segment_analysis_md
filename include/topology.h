// readfile.h
#ifndef _TOPOLOGY_H_
#define _TOPOLOGY_H_

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


/* ---------------------------------------------------------
   function declarations
   --------------------------------------------------------- */
/* void ReadMolecularContent2 (char *fname); */
/* void ReadLines2 (Chain *chain_segment); */
/* void topo_build(Chain *chain); */
void topo_build_intra(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],double cutoff);
void topo_build_inter(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS],double cutoff);

void topo_sort_dimer_contacts_initially(Chain *chain,int chain_num,FILE *fp); // dimer
void topo_sort_dimer_contacts_initially(Chain *chain,Chain *chain_ref, \
                                        int cid, FILE *fp,std::vector< std::vector <int> > contacts);

/* std::vector< std::vector <int> > pf_array; */

// not using CUTOFF ?
// overloading.
/* int topo_contacts_persisting(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS],double cutoff); */
int topo_contacts_persisting(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS]); // inter comparison
int topo_contacts_persisting(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],FILE *fp); // intra comparison
int topo_contacts_persisting(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS],double cutoff); // with cutoff
/* int topo_contacts_persisting(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS]);  // intra comparison */


// residue and its contacts
int topo_contacts_persisting_by_residue(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],FILE *fp);



void topo_clear_map(Contact (*map)[MAX_CONTACTS]);
int topo_count_map(Contact (*map)[MAX_CONTACTS]);

#endif
