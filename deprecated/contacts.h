#ifndef _CONTACTS_H_
#define _CONTACTS_H_
// contacts.h

/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "chain.h"

/* #include <stdio.h> */
/* #include <string> */
/* #include <stdlib.h> */
/* #include <iostream> */
/* #include <vector> */

// my headers
/* #include <md.h> */
#include <debug.h>
/* #include "dcdio.h" */


/* ---------------------------------------------------------
   definitions
   --------------------------------------------------------- */

/* ---------------------------------------------------------
   declarations
   --------------------------------------------------------- */
void test_contacts();
void build_map();
/* void build_map_intra(Chain *chain,Contact (*map)[MAX_CONTACTS],float cutoff); */
/* void build_map_intra(Chain *chain,Contact (*map)[MAX_CONTACTS],float cutoff); */
void build_map_intra(Chain *chain);

void get_overlapping_members();

#endif
