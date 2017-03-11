// readfile.h
#ifndef _TOPOLOGY_CHARMM_H_
#define _TOPOLOGY_CHARMM_H_

// libraries:
/* #include <stdio.h> */
/* #include <stdlib.h> */
/* #include <string> */

// Definitions:
/* #define BUFFERSIZE 900 */


/* ---------------------------------------------------------
   Class declarations
   --------------------------------------------------------- */
/* #ifndef CHAIN_H */
/* #define CHAIN_H */
/* #include "chain.h" */
/* #endif */


/* ---------------------------------------------------------
   function declarations
   --------------------------------------------------------- */
/* void topo_build_intra(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],double cutoff); */
/* void topo_build_inter(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS],double cutoff); */
void top_charmm_write_general(Chain *chain,int max_num_chains,FILE *fp);

#endif
