// readfile.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h> // strtod?, stod
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
extern int endoffile;

/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "readfile.h"
#include "chain.h"
#include "md.h"

// Vector pos;

// #define MOLECULE_INITIAL_SIZE 700
// Chain chain;


/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
void topo_sort_dimer_contacts_initially(Chain *chain, int cid, FILE *fp) {

    // WORKING, uses declaration inside of Class Chain for counting

    debug("hello from topo_dimer!!!\n");
    debug("%d %d \n",chain[cid].num_atoms_ca,cid);
    debug("northern:%d \n",chain[cid].LonN);

    // printf("NESW: %d %d %d %d\n",chain[cid].LonN, \
    //        chain[cid].LatE,chain[cid].LonS,chain[cid].LatW);

    // Doesn't necessarily work, because it wasn't tallied for chain_ref.
    // printf("WSEN %d< %d, %d> %d^\n", \
    //        chain[cid].LatW_total2x, \
    //        chain[cid].LonS_total2x, \
    //        chain[cid].LatE_total2x, \
    //        chain[cid].LonN_total2x);

    int alphaW, alphaS, alphaE, alphaN;
    alphaW = alphaS = alphaE = alphaN = 0;

    if (chain[cid].LatW != -1){
        alphaW = chain[cid].countcontacts(chain[cid].contactsLatW);
    }
    if (chain[cid].LonS != -1){
        alphaS = chain[cid].countcontacts(chain[cid].contactsLonS);
    }
    if (chain[cid].LatE != -1){
        alphaE = chain[cid].countcontacts(chain[cid].contactsLatE);
    }
    if (chain[cid].LonN != -1){
        alphaN = chain[cid].countcontacts(chain[cid].contactsLonN);
    }

    int beta;
    beta = chain[cid].LonN;

    int betaW, betaS, betaE, betaN;
    betaW = betaS = betaE = betaN = 0;

    // if (beta != -1) {
    if (chain[beta].LatW != -1){
        betaW = chain[beta].countcontacts(chain[beta].contactsLatW);
    }
    if (chain[beta].LonS != -1){
        betaS = chain[beta].countcontacts(chain[beta].contactsLonS);
    }
    if (chain[beta].LatE != -1){
        betaE = chain[beta].countcontacts(chain[beta].contactsLatE);
    }
    if (chain[beta].LonN != -1){
        betaN = chain[beta].countcontacts(chain[beta].contactsLonN);
    }

    // }
    // printf("WSEN %d< %d, %d> %d^\n",alphaW,alphaS,alphaE,alphaN);

#ifndef NDEBUG
    printf("%3d -S<>N:        %3d\n",alphaS,alphaN);
    printf("%3d      :S<>N--  %3d\n",betaS,betaN);
    printf("W(alpha)--%3d    %3d--W(beta)\n",alphaW,betaW);
    printf("E(alpha)--%3d    %3d--E(beta)\n",alphaE,betaE);
#endif

    int western, eastern, internal, southern, northern, externals;
    // internal: the alpha-beta monomer interfacial set of contacts, intradimer
    western = alphaW + betaW;
    eastern = alphaE + betaE;
    southern = alphaS;
    northern = betaN;
    externals = western + eastern + southern + northern;

    if (alphaN == betaS) {
        internal = alphaN;
    } else {
        internal = -1;
    }

    // fprintf(fp,"%d %d %d %d %d\n",southern,internal,northern,eastern,western);
    fprintf(fp,"%3d %6d %6d %3d %6d %6d   %3d %3d %3d %3d %3d %3d\n", \
            cid,chain[cid].index,chain[cid].findex, \
            beta,chain[beta].index,chain[beta].findex, \
            southern,internal,northern,eastern,western,externals);

}
// 2nd
void topo_sort_dimer_contacts_initially(Chain *chain,Chain *chain_ref, \
                                        int cid, FILE *fp, \
                                        std::vector< std::vector<int> > contacts){

// #ifndef NDEBUG
//     // check the 2D array: pf_array.
//     printf("insidetopo: 2d_array:\n");
//     for (int i=0; i<contacts.size(); i++) {
//         printf("%d:\t",i);
//         for (int j=0; j<contacts[i].size(); j++) {
//             printf("%3d ",contacts[i][j]);
//         }
//         printf("\n");
//     }
// #endif

    // printf("The size of contacts is double array: %d\n",contacts.size());
    // return;

    int alphaW, alphaS, alphaE, alphaN;
    alphaW = alphaS = alphaE = alphaN = 0;

    if (chain_ref[cid].LonN != -1){
        alphaN = contacts[cid][0];
    }
    if (chain_ref[cid].LatE != -1){
        alphaE = contacts[cid][1];
    }
    if (chain_ref[cid].LonS != -1){
        alphaS = contacts[cid][2];
    }
    if (chain_ref[cid].LatW != -1){
        alphaW = contacts[cid][3];
    }

    int beta;
    beta = chain_ref[cid].LonN;

    int betaW, betaS, betaE, betaN;
    betaW = betaS = betaE = betaN = 0;

    if (chain_ref[beta].LonN != -1){
        betaN = contacts[beta][0];
    }
    if (chain_ref[beta].LatE != -1){
        betaE = contacts[beta][1];
    }
    if (chain_ref[beta].LonS != -1){
        betaS = contacts[beta][2];
    }
    if (chain_ref[beta].LatW != -1){
        betaW = contacts[beta][3];
    }


#ifndef NDEBUG
    printf("%3d -S<>N:        %3d\n",alphaS,alphaN);
    printf("%3d      :S<>N--  %3d\n",betaS,betaN);
    printf("W(alpha)--%3d    %3d--W(beta)\n",alphaW,betaW);
    printf("E(alpha)--%3d    %3d--E(beta)\n",alphaE,betaE);
#endif

    int western, eastern, internal, southern, northern, externals;
    // internal: the alpha-beta monomer interfacial set of contacts, intradimer
    western = alphaW + betaW;
    eastern = alphaE + betaE;
    southern = alphaS;
    northern = betaN;
    externals = western + eastern + southern + northern;

    if (alphaN == betaS) {
        internal = alphaN;
    } else {
        internal = -1;
    }

    fprintf(fp,"%3d %6d %6d %3d %6d %6d   %3d %3d %3d %3d %3d %3d\n", \
            cid,chain[cid].index,chain[cid].findex, \
            beta,chain[beta].index,chain[beta].findex, \
            southern,internal,northern,eastern,western,externals);
}

// void topo_sort_dimer_contacts_later(FILE *fp,int) {
//     // later
// }

void topo_build_intra(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],double cutoff) {

    // printf("building topology [intra] for: %d at %f cutoff\n",chain[chain_num].chainid,cutoff);

    double dist = 0.0;
    int count = 0;

    // using 0,3 indices. 1,4, 2,5, ... 430,433
    for(int i=0; i<chain[chain_num].num_atoms_ca; i++){
        for(int j=i+3; j<chain[chain_num].num_atoms_ca; j++){
            dist = distance(chain[chain_num].pos[i],chain[chain_num].pos[j]);
            if(dist<cutoff){
                count ++;
                // printf("dist: %f\n",dist);
                // in resid i.
                for(int k=0; k<MAX_CONTACTS; k++){
                    if(map[i][k].cresid == -1){
                        map[i][k].cresid = j;
                        map[i][k].index = chain[chain_num].indices[i];
                        map[i][k].cindex = chain[chain_num].indices[j];
                        // double eh; // energy scaling (well depth)
                        map[i][k].distance = dist;
                        break;
                    } // if
                } // for 3
            } // if
        } // for 2
    } // for 1
    printf("total contacts [intra]: %d\n",count);
    // map[chain_num].total = count;
    return;
}

void topo_build_inter(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS],double cutoff) {
    // num1,num2: chains being compared.

    if(num2 == -1){
        debug("returning from %d -- %d\n",num1,num2);
        return;
    }

    // printf("building topology [inter] for: %d in contact with %d at %f cutoff\n",\
    //        chain[num1].chainid,chain[num2].chainid,cutoff);

    double dist = 0.0;
    int count = 0;



    // for(int i=0; i<chain[num1].num_atoms_ca; i++){
    //     printf("atom: %d\n",i);
    // }


    // using 0,3 indices. 1,4, 2,5, ... 430,433
    for(int i=0; i<chain[num1].num_atoms_ca; i++){
        for(int j=0; j<chain[num2].num_atoms_ca; j++){
            dist = distance(chain[num1].pos[i],chain[num2].pos[j]);
            if(dist<cutoff){
                count ++;
                // printf("dist: %f\n",dist);
                // in resid i.
                for(int k=0; k<MAX_CONTACTS; k++){
                    if(map[i][k].cresid == -1){
                        map[i][k].cresid = j;
                        map[i][k].index = chain[num1].indices[i];
                        map[i][k].cindex = chain[num1].indices[j];
                        // double eh; // energy scaling (well depth)
                        map[i][k].distance = dist;
                        break;
                    } // if
                } // for 3
            } // if
        } // for 2
    } // for 1
    debug("total contacts [inter]: %d\n",count);
    // // map[chain_num].total = count;

    return;
}

// overloaded
// int topo_contacts_persisting(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],double cutoff){
// int topo_contacts_persisting(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS]){
int topo_contacts_persisting(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],FILE *fp){
    // intra case
    // hence, the only chain_num

    std::vector<int> vecint_count (chain[chain_num].num_atoms_ca,0); // four ints with value 100

    double dist = 0.0;
    double dist_orig = 0.0;
    int count = 0; // local (per residue)
    int total = 0; // total (per chain)

    // using 0,3 indices. 1,4, 2,5, ... 430,433
    // for(int i=0; i<chain[chain_num].num_atoms_ca; i++){
        // for(int j=i+3; j<chain[chain_num].num_atoms_ca; j++){
    for(int i=0; i<chain[chain_num].num_atoms_ca; i++){

        count = 0;
        for(int k=0; k<MAX_CONTACTS; k++){

            if(map[i][k].cresid == -1){
                break;
            } else {

                dist = distance(chain[chain_num].pos[i],chain[chain_num].pos[map[i][k].cresid]);
                // printf("the distance was: %f, it is now: %f\n",map[i][k].distance,dist);
                dist_orig = map[i][k].distance;

                if(dist < dist_orig + 2.0  || dist < 8.0 ){
                    vecint_count[i] ++;
                    vecint_count[map[i][k].cresid] ++;
                    count ++;
                    total ++;
                }
            } // if
        } // for k
    } // for i

    // NEW
    // for(int x=0; x<chain[chain_num].num_atoms_ca; x++){
    //     // debug("%d\n",vecint_count[x]);
    //     fprintf(fp,"%d ",vecint_count[x]);
    // }

    int sum_of_elems = 0;
    for(std::vector<int>::iterator it = vecint_count.begin(); it != vecint_count.end(); ++it){
        sum_of_elems += *it;
        // fprintf(fp,"%d ",vecint_count[x]);
        fprintf(fp,"%d ",*it);
    }

    printf("total contacts that persisted [intra]: %d   %d  (using 2x)\n",total,sum_of_elems);
    fprintf(fp,"\n# remaining in frame(above): %d   %d (<-- used double counted)\n",total,sum_of_elems);
    return count;
}
int topo_contacts_persisting(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS]){
    // num1,num2: chains being compared.
    // "original distance" is created from using the chain_ref Contact map.

    if(num2 == -1){
        debug("returning from %d -- %d\n",num1,num2);
        return 0;
    }
    // printf("checking topology [inter] for: %d in contact with %d at %f cutoff\n",\
    //        chain[num1].chainid,chain[num2].chainid,cutoff);

    double dist = 0.0;
    double dist_orig = 0.0;
    int count = 0;

    for(int i=0; i<chain[num1].num_atoms_ca; i++){
        for(int k=0; k<MAX_CONTACTS; k++){

            if(map[i][k].cresid == -1){
                break;
            } else {
                dist = distance(chain[num1].pos[i],chain[num2].pos[map[i][k].cresid]);
                // printf("the distance was: %f, it is now: %f\n",map[i][k].distance,dist);
                dist_orig = map[i][k].distance;

                // PERSISTENCE CRITERION!
                // 8.0 to 20.0 hard definition
                // if(dist < cutoff){
                //     count ++;
                // }
                if(dist < dist_orig + 2.0  || dist < 8.0 ){
                    count ++;
                }
            }
        } // for 2
    } // for 1
    debug("total contacts that persisted: %d\n",count);
    return count;
}
int topo_contacts_persisting(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS],double cutoff)
{
    // num1,num2: chains being compared.
    // "original distance" is created from using the chain_ref Contact map.

    if(num2 == -1){
        debug("returning from %d -- %d\n",num1,num2);
        return 0;
    }
    // printf("checking topology [inter] for: %d in contact with %d at %f cutoff\n",\
    //        chain[num1].chainid,chain[num2].chainid,cutoff);

    double dist = 0.0;
    double dist_orig = 0.0;
    int count = 0;

    for(int i=0; i<chain[num1].num_atoms_ca; i++){
        for(int k=0; k<MAX_CONTACTS; k++){

            if(map[i][k].cresid == -1){
                break;
            } else {
                dist = distance(chain[num1].pos[i],chain[num2].pos[map[i][k].cresid]);
                // printf("the distance was: %f, it is now: %f\n",map[i][k].distance,dist);
                dist_orig = map[i][k].distance;

                // PERSISTENCE CRITERION!
                // 8.0 to 20.0 hard definition
                if(dist < cutoff){
                    count ++;
                }
                // if(dist < dist_orig + 2.0  || dist < 8.0 ){
                //     count ++;
                // }
            }
        } // for 2
    } // for 1
    debug("total contacts that persisted: %d\n",count);
    return count;
}

int topo_contacts_persisting_by_residue(Chain *chain,int chain_num,Contact (*map)[MAX_CONTACTS],FILE *fp){
    // intra case
    // hence, the only chain_num

    std::vector<int> vecint_count (chain[chain_num].num_atoms_ca,0); // four ints with value 100

    double dist = 0.0;
    double dist_orig = 0.0;
    int count = 0; // local (per residue)
    int total = 0; // total (per chain)

    // using 0,3 indices. 1,4, 2,5, ... 430,433
    // for(int i=0; i<chain[chain_num].num_atoms_ca; i++){
        // for(int j=i+3; j<chain[chain_num].num_atoms_ca; j++){
    for(int i=0; i<chain[chain_num].num_atoms_ca; i++){

        // FPRINT RES HERE
        // printf("%d\n",i);
        fprintf(fp,"%3d ",i);

        count = 0;
        for(int k=0; k<MAX_CONTACTS; k++){

            if(map[i][k].cresid == -1){

                fprintf(fp,"%3d ",-1);
                // break;
            } else {

                dist = distance(chain[chain_num].pos[i],chain[chain_num].pos[map[i][k].cresid]);
                // printf("the distance was: %f, it is now: %f\n",map[i][k].distance,dist);
                dist_orig = map[i][k].distance;

                if(dist < dist_orig + 2.0  || dist < 8.0 ){
                    vecint_count[i] ++;
                    vecint_count[map[i][k].cresid] ++;
                    // vecint_count[map[]]
                    count ++;
                    total ++;

                    // FPRINT RES CONTACT HERE
                    // printf("\t%d\n",map[i][k].cresid);
                    fprintf(fp,"%3d ",map[i][k].cresid);

                }
                else
                {
                    fprintf(fp,"%3d ",-1);
                }
            } // if
        } // for k

        fprintf(fp,"\n");
    } // for i

    // NEW
    // for(int x=0; x<chain[chain_num].num_atoms_ca; x++){
    //     // debug("%d\n",vecint_count[x]);
    //     fprintf(fp,"%d ",vecint_count[x]);
    // }

    // int sum_of_elems = 0;
    // for(std::vector<int>::iterator it = vecint_count.begin(); it != vecint_count.end(); ++it){
    //     sum_of_elems += *it;
    //     // fprintf(fp,"%d ",vecint_count[x]);
    //     // fprintf(fp,"%d ",*it);
    //     // printf("%d\n",*it);
    //     // fprintf("%d\n",it - vecint_count.begin());
    //     printf("%d\n",it - vecint_count.begin());

    //     // // 20 positions + 1
    //     // for(int res; res < 21; res++)
    //     // {

    //     // }
    // }
    // printf("total contacts that persisted [intra]: %d   %d  (using 2x)\n",total,sum_of_elems);
    // fprintf(fp,"\n# remaining in frame(above): %d   %d (<-- used double counted)\n",total,sum_of_elems);

    return count;
}


int topo_contactmap_persisting(Chain *chain,int num1,int num2,Contact (*map)[MAX_CONTACTS]){
    // num1,num2: chains being compared.
    // "original distance" is created from using the chain_ref Contact map.

    // map of remaining contacts.
    Contact rem[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];

    if(num2 == -1){
        debug("returning from %d -- %d\n",num1,num2);
        return 0;
    }
    // printf("checking topology [inter] for: %d in contact with %d at %f cutoff\n",\
    //        chain[num1].chainid,chain[num2].chainid,cutoff);

    double dist = 0.0;
    double dist_orig = 0.0;
    int count = 0;

    for(int i=0; i<chain[num1].num_atoms_ca; i++){
        for(int k=0; k<MAX_CONTACTS; k++){

            if(map[i][k].cresid == -1){
                break;
            } else {
                dist = distance(chain[num1].pos[i],chain[num2].pos[map[i][k].cresid]);
                // printf("the distance was: %f, it is now: %f\n",map[i][k].distance,dist);
                dist_orig = map[i][k].distance;

                // PERSISTENCE CRITERION!
                // 8.0 to 20.0 hard definition
                // if(dist < cutoff){
                //     count ++;
                // }
                if(dist < dist_orig + 2.0  || dist < 8.0 ){
                    count ++;
                }
            }
        } // for 2
    } // for 1
    debug("total contacts that persisted: %d\n",count);
    return count;
}

int topo_count_map(Contact (*map)[MAX_CONTACTS]) {
    // num1,num2: chains being compared.

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
        debug("total contacts in this map are: %d\n",count);
    } else {
        debug("empty map.\n");
    }

    // // map[chain_num].total = count;

    return count;
}

void topo_clear_map(Contact (*map)[MAX_CONTACTS]){

    for(int i=0; i<MOLECULE_INITIAL_SIZE; i++){
        for(int k=0; k<MAX_CONTACTS; k++){
            map[i][k].cresid = -1;
            map[i][k].index = -1;
            map[i][k].cindex = -1;
            map[i][k].total = -1;
            map[i][k].eh = 0.0;
            map[i][k].distance = 0.0;
            // map[i][k].~Contact(); ** adds contacts, round to round.
        }
    }
}
