// readfile.h
#ifndef _READFILE_H_
#define _READFILE_H_

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
#include "chain.h"
/* #endif */


/* ---------------------------------------------------------
   function declarations
   --------------------------------------------------------- */
void ReadMolecularContent (char *fname);
void ReadLines (Chain *chain_segment);
void ReadEveryLine (Chain *chain);
#endif
