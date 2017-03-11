// chi.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h> // strtod?, stod
#include <math.h>

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



/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
double get_chi_global(Chain *chain_ref,Chain *chain_later,int L1,int L2){

    debug("getting chi..\n");



    // int N = 0;
    double N,prefactor,theta,chi,distO,distN,diff;
    N = prefactor = theta = chi = distO = distN = diff = 0.0;

    // N = (double) chain_ref[0].num_atoms_ca;
    N = abs((double) L1 - (double) L2) + 1.0;
    debug("N: %2.10f\n",N);
    prefactor = 2 / (N*N - 5 * N + 6);

    debug("prefactor: %f\n",prefactor);


    // 2.0 should be def_param defined.
    //             if(sum <= R_limit_new) diff ++;

    for(int i=L1; i<=(L2-3); i++)
    {

        for(int j=(i+3); j<=L2; j++)
        {

            // distO: original, distN: Now.
            distO = distance(chain_ref[0].pos[i],chain_ref[0].pos[j]);
            distN = distance(chain_later[0].pos[i],chain_later[0].pos[j]);
            diff = fabs(distN - distO);

            // 2.0 possibly the default for epsilon
            if(diff <= 2.0){
                theta += 1.0;
            }
        }
    }
    chi = 1 - prefactor * theta;
    debug("chi: %f\n",chi);

    return chi;
}

void get_chi_by_residue(Chain *chain_ref,Chain *chain_later,int L1,int L2,FILE *fp){

    debug("getting chi_by_residue.. %d-%d\n",L1,L2);
    // double N,prefactor,theta,chi,distO,distN,diff;
    // N = prefactor = theta = chi = distO = distN = diff = 0.0;

    int total,diff;
    double distO,distN;
    double sum,chi;
    total = 0;
    diff = 0;
    sum  = 0.0;
    chi = 0.0;

    /* ------------------- get the Chi value for the stretch between P1 and P2 ------------------- */
    for(int i=L1;i<=L2;i++)
    {
        // for(i = P1 ; i <= P2 ; i ++){
        //for(i = 199; i <= 296; i++){
        total = 0;
        diff = 0;
        sum  = 0.0;
        chi = 0.0;

        for(int j=L1;j<=L2;j++)
        {
            //for(j = 199; j <= 296; j++){

            if(abs(i-j)>=2)
            {
                // sum = fabs(dist[i][j] - dist0[i][j]);
                // if(sum <= R_limit_new) diff ++;
                // total ++;
                //fprintf(out,"%5d %5d %10d\n", i, j, diff);

                // distO: original, distN: Now.
                distO = distance(chain_ref[0].pos[i],chain_ref[0].pos[j]);
                distN = distance(chain_later[0].pos[i],chain_later[0].pos[j]);
                sum = fabs(distN - distO);

                if(sum <= 2.0)
                {
                    diff ++;
                }
                total ++;
            }
        }
        chi = 1.0 - ((double) diff)/((double) total);

        /* commands for individual writing */
        fprintf(fp,"%4.3f ",chi);

        // sprintf(filename_write, "Chi_residue_%d.dat", i);
        // out_chi = fopen(filename_write, "a");
        // fprintf(out_chi,"%10.3f \n", Chi);
        // fclose(out_chi);

    } // FOR

    fprintf(fp,"\n",chi);

}
