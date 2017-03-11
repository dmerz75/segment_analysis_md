// main.cpp

// headers C
extern "C" {
// your functions here for the header

}

// my headers
#include "debug.h"
#include "readfile.h"
#include "chain.h"
#include "md.h"
// #include "contacts.h"
#include "topology.h"
#include "topology_charmm.h"
#include "curvature.h"
#include "chi.h"
#include "dcd.h"

// headers C++
#include <stdlib.h>
#include <stdio.h> // printf
#include <string.h> // strcpy, memcpy
#include <new> // delete
#include <ctype.h> // getopt - stuff
#include <unistd.h> // getopt - stuff
#include <vector>
#include <algorithm> // bool & sort
#include <iostream>
#include <iomanip> // setw
#include <fstream> //
#include <map> // map

// boost
// #include "boost/multi_array.hpp"

// VMD
// #ifndef _DCDIO_H_
// #define _DCDIO_H_
// #include "dcdio.h"
// #endif

// namespaces

int main(int argc, char *argv[]) {

    debug(">> debugging is turned on!\n");
    debug("Analysis by Segment, Beginning with %d arguments.\n-->\n",argc);
    printf("Welcome to the segmental/chain analyzer!\n");
    printf("chain allocation size is: %d\n",MOLECULE_INITIAL_SIZE);
    printf("chain allocation size is(mt): %d\n",MOLECULE_INITIAL_SIZE_MT);

    // START HERE.
    int num_chains=0, chains_ignore=0, chains_to_use=0;


#ifdef DCDREAD
    int frame_position;
    frame_position = 0;
    int start,stop,step;

    if (argc < 8) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference-PDB>" \
                  << " <Filename-timelater-DCD>" \
                  << " <num_chains>" \
                  << " <chains_ignore>" \
                  << " <start>" \
                  << " <stop>" \
                  << " <step>" \
                  <<std::endl;
        exit(1);
    }
    // run_segment  ref   timelater maxchains ignorechains
    // argc 1        2          3        4           5
    // argv 0        1          2        3           4
    num_chains = atoi(argv[3]); // argv = 3, chains are 0,1,2
    chains_ignore = atoi(argv[4]);
    chains_to_use = num_chains - chains_ignore;

    // start,stop,step
    start = atoi(argv[5]);
    stop = atoi(argv[6]);
    step = atoi(argv[7]);


    // REFERENCE: chain_pdb (PDB-Big)
    Chain *chain_pdb;
    try {
        chain_pdb = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }
    for ( int i = 0; i < num_chains; i ++ ) {
        debug("__building chain --> %d\n",i);
        strcpy(chain_pdb[i].filename,argv[1]);
        chain_pdb[i].chainid = i;
        if ( i == 0 ) {
            // debug("__building first chain (%d)\n", i);
            chain_pdb[i].index = 0;
            chain_pdb[i].file_line_begin = 0;
        } else {
            // debug("__building more chains --> %d\n",i);
            chain_pdb[i].index = chain_pdb[i-1].findex + 1;
            chain_pdb[i].file_line_begin = chain_pdb[i-1].file_line_end; // or +1?
        }
        // Read PDB File.
        ReadLines(&chain_pdb[i]);
        chain_pdb[i].findex = chain_pdb[i].index + chain_pdb[i].num_atoms_ca - 1;
        chain_pdb[i].assign_indices();

        if (chain_pdb[i].num_atoms_ca == 0) {
            num_chains = i;
            debug("Reached end of chain.\n");
            break;
        }

#ifdef INFO
        if (chain_pdb[i].num_atoms_ca > 0) {chain_pdb[i].print_prop();}
#endif // INFO
        // if (chain_pdb[i].num_atoms_ca > 0) {chain_pdb[i].print_prop();}
    }
    // num_chains = i + 1;
    std::cout << " <<<<<<<< chain_pdb has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;
    // End of REFERENCE chain_pdb.
    // exit(0);



    // REFERENCE: chain_ref (PDB-Small)
    Chain *chain_ref;
    try {
        chain_ref = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }
    for ( int i = 0; i < num_chains; i ++ ) {
        // Undefine.
        debug("Redefining molecule size: (%d)\n",chain_pdb[i].num_atoms_ca);

        // RIGHT SIZE !!!  <-- Read again. (know num_atoms_ca).
#undef MOLECULE_INITIAL_SIZE
#define MOLECULE_INITIAL_SIZE chain_pdb[i].num_atoms_ca

        debug("__building chain --> %d\n",i);
        strcpy(chain_ref[i].filename,argv[1]);
        chain_ref[i].chainid = i;
        if ( i == 0 ) {
            debug("__building first chain (%d)\n", i);
            chain_ref[i].index = 0;
            chain_ref[i].file_line_begin = 0;
        } else {
            debug("__building more chains --> %d\n",i);
            chain_ref[i].index = chain_ref[i-1].findex + 1;
            chain_ref[i].file_line_begin = chain_ref[i-1].file_line_end; // or +1?
        }
        // Read PDB File.
        ReadLines(&chain_ref[i]);
        chain_ref[i].findex = chain_ref[i].index + chain_ref[i].num_atoms_ca - 1;
        chain_ref[i].assign_indices();
#ifdef INFO
        if (chain_ref[i].num_atoms_ca > 0) {chain_ref[i].print_prop();}
#endif // INFO
    }
    std::cout << " <<<<<<<< chain_ref has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;
    // End of REFERENCE chain_ref



    /* ---------------------------------------------------------
       Deleted chain_pdb.
       --------------------------------------------------------- */
    debug("removing chain_pdb.\n");
    delete [] chain_pdb;
    /* ---------------------------------------------------------
       Deleted chain_pdb.
       --------------------------------------------------------- */


    // REFERENCE: chain_0 (DCD-0)
    Chain *chain_0;
    try {
        chain_0 = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }
    for ( int i = 0; i < num_chains; i ++ ) {
        // Undefine.
        // debug("Redefining molecule size: (%d)\n",chain_ref[i].num_atoms_ca);
#undef MOLECULE_INITIAL_SIZE
#define MOLECULE_INITIAL_SIZE chain_ref[i].num_atoms_ca

        debug("__building chain --> %d\n",i);
        strcpy(chain_0[i].filename,argv[1]);
        chain_0[i].chainid = i;
        if ( i == 0 ) {
            debug("__building first chain (%d)\n", i);
            chain_0[i].index = 0;
            chain_0[i].file_line_begin = 0;
        } else {
            debug("__building more chains --> %d\n",i);
            chain_0[i].index = chain_0[i-1].findex + 1;
            chain_0[i].file_line_begin = chain_0[i-1].file_line_end; // or +1?
        }
        // Read PDB File.
        ReadLines(&chain_0[i]);
        chain_0[i].findex = chain_0[i].index + chain_0[i].num_atoms_ca - 1;
        chain_0[i].assign_indices();
#ifdef INFO
        if (chain_0[i].num_atoms_ca > 0) {chain_0[i].print_prop();}
#endif // INFO
    }
    std::cout << " <<<<<<<< chain_0 has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;
    // End of REFERENCE chain_0


    // REFERENCE: chain_later
    Chain *chain_later;
    try {
        chain_later = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }
    for ( int i = 0; i < num_chains; i ++ ) {
        // Undefine.
        debug("Redefining molecule size: (%d)\n",chain_ref[i].num_atoms_ca);

        // RIGHT SIZE !!!  <-- Read again. (know num_atoms_ca).
#undef MOLECULE_INITIAL_SIZE
#define MOLECULE_INITIAL_SIZE chain_ref[i].num_atoms_ca

        debug("__building chain --> %d\n",i);
        strcpy(chain_later[i].filename,argv[1]);
        chain_later[i].chainid = i;
        if ( i == 0 ) {
            debug("__building first chain (%d)\n", i);
            chain_later[i].index = 0;
            chain_later[i].file_line_begin = 0;
        } else {
            debug("__building more chains --> %d\n",i);
            chain_later[i].index = chain_later[i-1].findex + 1;
            chain_later[i].file_line_begin = chain_later[i-1].file_line_end; // or +1?
        }
        // Read PDB File.
        ReadLines(&chain_later[i]);
        chain_later[i].findex = chain_later[i].index + chain_later[i].num_atoms_ca - 1;
        chain_later[i].assign_indices();
#ifdef INFO
        if (chain_later[i].num_atoms_ca > 0) {chain_later[i].print_prop();}
#endif // INFO
    }
    std::cout << " <<<<<<<< chain_later has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;
    // End of REFERENCE chain_later

    /* ---------------------------------------------------------
      Built chain_ref.
      Built chain_0.
      Built chain_later.
      --------------------------------------------------------- */


    /* ---------------------------------------------------------
       write dcd
       --------------------------------------------------------- */
#ifdef DCD_WRITE_B

    std::string str_dcd_read(argv[2]);
    std::size_t found;
    found=str_dcd_read.find(".dcd",0);

    std::string str_dcd_name = str_dcd_read.substr(0,found);
    std::string str_dcd_write = str_dcd_name + "_subset.dcd";
    char *fn_dcd_write = (char*) str_dcd_write.c_str();

    printf("\nreading dcd >>>  <%s>\n",str_dcd_read.c_str());
    printf("using name for dcd >>>  <%s>\n",str_dcd_name.c_str());
    printf("writing dcd >>>  <%s>\n",fn_dcd_write);

    // Write DCD
    molfile_timestep_t timestep_w;
    void *vw;
    dcdhandle *dcdw;
    int natoms_w = 0;

    dcdw = (dcdhandle *)vw;

    // get atom total for writing.
    for(int i=0; i<num_chains; i++){
        natoms_w += chain_ref[i].num_atoms_ca;
    }
    printf("for dcd writing >>>  <%d> atoms expected.\n",natoms_w);


    vw = open_dcd_write(fn_dcd_write,"dcd",natoms_w);
    if (!vw) {
        fprintf(stderr, "main) open_dcd_write failed for file %s\n", *fn_dcd_write);
        return 1;
    } else {
        printf("opened <%s> successfully!!\n\n",fn_dcd_write);
    }

    timestep_w.coords = (float *)malloc(3*sizeof(float)*natoms_w);

    // dcd = (dcdhandle *)v;
    // sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    // totalMB += sizeMB;
    // printf("main) file: %s\n", *argv);
    // printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);
    // timestep.coords = (float *)malloc(3*sizeof(float)*natoms);

#endif // DCD_WRITE_B



    /* ---------------------------------------------------------
       LatLon Neighbors
       --------------------------------------------------------- */
// #if defined(LATLON) || defined(LONONLY_PROTO)
// #endif // #if defined(LATLON) || defined(LONONLY_PROTO)

#ifdef LATLON
    std::cout << "getting lateral and longitudinal neighbors " << chains_to_use << std::endl;

    Vector centroid_front;
    Vector centroid_back;
    Vector mt_axis;
    Vector mt_axis_norm;

    // ComputeCentroid()
    for ( int i=0; i < num_chains; i++) {
        chain_ref[i].ComputeCentroid();
    }

    // get front 10%, back 10%
    int front10, back10 = 0;
    front10 = (int) (chains_to_use * 0.10);
    back10 = (int) (chains_to_use * 0.90);

    printf("front10 and back10: %d %d\n",front10,back10);

    // centroid_front = get_centroid_from_centroids(chain_ref, 10, 40);
    centroid_front = get_centroid_from_centroids(chain_ref, 0, front10);
    // centroid_back = get_centroid_from_centroids(chain_ref, chains_to_use - 46, chains_to_use-15);
    centroid_back = get_centroid_from_centroids(chain_ref,back10, chains_to_use);
    // exit(0);

    // print centroids (2)
    get_vector(centroid_front, centroid_back, &mt_axis);
    normalize(mt_axis, &mt_axis_norm);
    debug("Centroid_front: %f %f %f\n",centroid_front.x,centroid_front.y,centroid_front.z);
    debug("Centroid_back: %f %f %f\n",centroid_back.x,centroid_back.y,centroid_back.z);
    debug("mt_axis: %f %f %f\n",mt_axis.x,mt_axis.y,mt_axis.z);
    debug("mt_axis_norm: %f %f %f\n",mt_axis_norm.x,mt_axis_norm.y,mt_axis_norm.z);

    printf(".. computing a general set of nearest neighbors .......\n");
    compute_nearest_neighbors(chain_ref,chains_to_use);

    printf(".. determining the 4 lateral and longitudinal neighbors .......\n");
    determine_latlon_neighbors(chain_ref,chains_to_use,mt_axis_norm);
#endif // LATCON

#ifdef MTPF_B
    printf("MTPF_B: processing..\n");

    std::vector<int> pf_south;
    std::vector<int>::iterator it; // monomers with no south
    std::vector<int>::iterator it2; // northern neighbors for any 1 monomer with no south.
    // std::vector<int>::iterator m1;
    // std::vector<int>::iterator m2;
    pf_south = get_monomers_with_no_southern_neighbor(chain_ref,chains_to_use);
    printf("size: %d\n",pf_south.size()); // usually 13.
    // printf("(main) protofilaments identified:\n");
    // for (it=pf_south.begin(); it<pf_south.end(); it++) {
    //     printf("%d ",*it);
    // }
    // std::cout << ' ' << std::endl; // has the implicit \n


    std::vector< std::vector <int> > pf_array;
    std::vector<int> pf_monomers;

    for (it=pf_south.begin(); it<pf_south.end(); it++) {
        // printf("%d ",*it);
        printf("%d\n",*it);
        pf_monomers = get_protofilament(chain_ref,*it);

        printf("(main) protofilament's neighbors identified:\n");
        for (it2=pf_monomers.begin(); it2<pf_monomers.end(); it2++) {
            printf("%d ",*it2);
        }
        std::cout << ' ' << std::endl; // has the implicit \n

        pf_array.push_back(pf_monomers); // ADD the ROW.
        pf_monomers.clear();
    }



    // check the 2D array: pf_array.
    printf("2d_array:\n");
    for (int i = 0; i<pf_array.size(); i++) {
        printf("%d:\t",i);
        for (int j = 0; j<pf_array[i].size(); j++) {
            printf("%3d ",pf_array[i][j]);
        }
        printf("\n");
    }




    printf("OPENING FILES:\n");
    std::ofstream outcontactfile[pf_array.size()*4];

    for(int i = 0; i < pf_array.size(); i++){
        std::ostringstream filename1;
        std::ostringstream filename2;
        std::ostringstream filename3;
        std::ostringstream filename4;

        filename1 << "pfcontacts_" << i << "_n.dat";
        filename2 << "pfcontacts_" << i << "_e.dat";
        filename3 << "pfcontacts_" << i << "_s.dat";
        filename4 << "pfcontacts_" << i << "_w.dat";

        std::cout << filename1.str().c_str() << std::endl;
        std::cout << filename2.str().c_str() << std::endl;
        std::cout << filename3.str().c_str() << std::endl;
        std::cout << filename4.str().c_str() << std::endl;

        outcontactfile[i*4].open(filename1.str().c_str());
        outcontactfile[i*4+1].open(filename2.str().c_str());
        outcontactfile[i*4+2].open(filename3.str().c_str());
        outcontactfile[i*4+3].open(filename4.str().c_str());
    }



    Contact map_ref[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
    int num_contacts;
    std::vector< std::vector <int> > max_contacts;
    std::vector<int> max_contacts_row; // N E S W
    // std::tuple<int,int,int> max_contacts();

    // ONLY going to run if System has contactsLonN, LatE, LonS, and LatW defined.
    // if(chains_to_use > 1){
    for (int i=0; i<chains_to_use; i++) {
        // for (int i=0; i<1; i++) {

        // N
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LonN,map_ref,8.0);
        memcpy(chain_ref[i].contactsLonN,map_ref, sizeof(chain_ref[i].contactsLonN));
        num_contacts = topo_count_map(map_ref);
        max_contacts_row.push_back(num_contacts);
        topo_clear_map(map_ref);
        // E
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LatE,map_ref,8.0);
        memcpy(chain_ref[i].contactsLatE,map_ref, sizeof(chain_ref[i].contactsLatE));
        num_contacts = topo_count_map(map_ref);
        max_contacts_row.push_back(num_contacts);
        topo_clear_map(map_ref);
        // S
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LonS,map_ref,8.0);
        memcpy(chain_ref[i].contactsLonS,map_ref, sizeof(chain_ref[i].contactsLonS));
        num_contacts = topo_count_map(map_ref);
        max_contacts_row.push_back(num_contacts);
        topo_clear_map(map_ref);
        // W
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LatW,map_ref,8.0);
        memcpy(chain_ref[i].contactsLatW,map_ref, sizeof(chain_ref[i].contactsLatW));
        num_contacts = topo_count_map(map_ref);
        max_contacts_row.push_back(num_contacts);
        topo_clear_map(map_ref);

        max_contacts.push_back(max_contacts_row);
        max_contacts_row.clear();
    }

#ifndef NDEBUG
    for ( int i = 0; i<max_contacts.size(); i++ ) {
        printf("contact: %d\n",i);
        for ( int d = 0; d < 4; d++ ) {
            printf("%d ",max_contacts[i][d]);
        }
        printf("\n");
    }
#endif


    for (int i = 0; i<pf_array.size(); i++) {

        outcontactfile[i*4+0] << "# contacts(n):\n" << "# chains: ";
        outcontactfile[i*4+1] << "# contacts(e):\n" << "# chains: ";
        outcontactfile[i*4+2] << "# contacts(s):\n" << "# chains: ";
        outcontactfile[i*4+3] << "# contacts(w):\n" << "# chains: ";

        // N: chainid(s)
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4] << pf_array[i][j] << " ";
        }
        outcontactfile[i*4] << "\n# max: ";
        // N: max_contacts
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4] << pf_array[i][j] << "|" << max_contacts[pf_array[i][j]][0] << " ";
        }
        outcontactfile[i*4] << "\n";


        // E: chainid(s)
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4+1] << pf_array[i][j] << " ";
        }
        outcontactfile[i*4+1] << "\n# max: ";
        // E: max_contacts
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4+1] << pf_array[i][j] << "|" << max_contacts[pf_array[i][j]][1] << " ";
        }
        outcontactfile[i*4+1] << "\n";


        // S: chainid(s)
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4+2] << pf_array[i][j] << " ";
        }
        outcontactfile[i*4+2] << "\n# max: ";
        // S: max_contacts
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4+2] << pf_array[i][j] << "|" << max_contacts[pf_array[i][j]][2] << " ";
        }
        outcontactfile[i*4+2] << "\n";


        // W: chainid(s)
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4+3] << pf_array[i][j] << " ";
        }
        outcontactfile[i*4+3] << "\n# max: ";
        // W: max_contacts
        for (int j = 0; j<pf_array[i].size(); j++ ) {
            outcontactfile[i*4+3] << pf_array[i][j] << "|" << max_contacts[pf_array[i][j]][3] << " ";
        }
        outcontactfile[i*4+3] << "\n";
    }


    printf("MTPF_B: complete.\n");
#endif // MTPF_B


#ifdef LONONLY_PROTO
    std::cout << "setting up protofilament analysis: angles and contacts, chains:" << chains_to_use << std::endl;

    // ComputeCentroid()
    for ( int i=0; i < num_chains; i++) {
        chain_ref[i].ComputeCentroid();
    }

    for(int i=0; i<chains_to_use; i++ ){

        if ( i != chains_to_use - 1) {
            chain_ref[i].LonN = i + 1;
        }

        if ( i != 0 ) {
            chain_ref[i].LonS = i - 1;
        }
    }
#endif // LONONLY_PROTO

#ifdef DIMERMAP_1A
    printf("beginning section dimermap_1A now ..\n");


    std::vector<int> dimermap_alphas;
    // for (int j=0; j<dimermap_alphas.size(); j++) {

    // build dimermap_alphas
    for (int i=0; i<chains_to_use; i++)
    {
        // alpha-tubulin
        if (chain_ref[i].num_atoms_ca >= 430)
        {
            dimermap_alphas.push_back(i);
        }
    }

    // alpha cid i f | beta cid i f | sinewe->12 | angles 15 | curv 19.
    // evaluation of: contacts(sinewe), angles(alat,blat,blon),curv(alatcur,blatcur,aloncur,bloncur)
    std::map <int , std::vector<int> > map_contacts; // 104x6
    std::map <int , std::vector<double> > map_angles; // 104x3
    std::map <int , std::vector<double> > map_curv; // 104x4

    std::vector<int> maprow_contacts; // nesw by monomer
    std::vector<double> maprow_angles; // by dimer
    std::vector<double> maprow_curv;

    // FILE
    FILE * fp_mt_analysis;
    fp_mt_analysis = fopen("mt_analysis.dat", "w+");

    fprintf(fp_mt_analysis,"# start_frame: %d  stop: %d  step: %d\n",start,stop,step);
    fprintf(fp_mt_analysis,"# 0-5   alpha index findex | beta index findex\n");
    fprintf(fp_mt_analysis,"# 6-11  SINEW: southern,internal,northern,eastern,western,external\n");
    fprintf(fp_mt_analysis,"# 12-14 angles: alpha lateral, beta lateral, beta longitudinal\n");
    // fprintf(fp_mt_analysis,"# 15-17|18-20 RofCurva: a_lat_Rcurv, b_lat_Rcurv, b_lon_Rcurv, a_lat_radius, b_lat_radius, b_lon_radius\n");
    fprintf(fp_mt_analysis,"# 15-17 RofCurva: a_lat_radius, b_lat_radius, b_lon_radius\n");

#endif // DIMERMAP_1A


#ifdef CENTROIDMOVEMENT_B

    // FILE
    FILE * fp_xyzpos;
    fp_xyzpos = fopen("centroid_xyzpos.dat", "a+");
    fprintf(fp_xyzpos,"# centroid_xyzpos.dat\n");
    // fprintf(fp_xyzpos,"\n");
    // fclose(fp_xyzpos);

    // stop
    // fprintf(fp_xyzpos,"# last_frame: %d\n",frame_position);
    // fprintf(fp_xyzpos,"# stop_frame: %d\n",stop);
    // fclose(fp_xyzpos);

#endif // CENTROIDMOVEMENT_E

#ifdef CHI_BEGIN

    // if (argc < 10) {
    if (argc < 8) {
        std::cout << "For chi: use 10 args." << std::endl;
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference-PDB>" \
                  << " <Filename-timelater-DCD>" \
                  << " <num_chains>" \
                  << " <chains_ignore>" \
                  << " <start>" \
                  << " <stop>" \
                  << " <step>" \
                  << " <limit-begin>" \
                  << " <limit-end>"
                  <<std::endl;
        exit(1);
    }
    // start,stop,step
    // start = atoi(argv[5]);
    // stop = atoi(argv[6]);
    // step = atoi(argv[7]);

    int limit_begin, limit_end;
    limit_begin = limit_end = -1;
    // limit_begin = atoi(argv[8]);
    // limit_end = atoi(argv[9]);

    limit_begin = 0;
    limit_end = chain_ref[0].num_atoms_ca-1;

    // FILE
    FILE * fp_chi;
    // char filename_chi[] = "chi_analysis.dat";
    // char filename_chi[40] = ("chi_analysis_%d_%d.dat",limit_begin,limit_end);
    // char filename_chi[] = ("chi_analysis_%s_%s.dat",argv[8],argv[9]);
    // std::string my_string1 = "a string %s %s", argv[8],argv[9];
    // std::string my_string1 = "a string %s %s",'hello','world';
    char filename_chi[sizeof "chi_analysis_100_100.dat"];
    sprintf(filename_chi,"chi_analysis_%d_%d.dat",limit_begin,limit_end);
    fp_chi = fopen(filename_chi,"w");

#endif // CHI BEGIN


#ifdef CURVATURE_B
    // usage
    if (argc != 8) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-PDB>" \
                  << " <Filename-DCD>" \
                  << " <num_chains>" \
                  << " <ignore_chains>" \
                  << " <start>" \
                  << " <stop>" \
                  << " <step>" \
                  << "\ncreates: distance_proto_centroids.dat (a+)" \
                  << "\ncreates: curvature_global.dat (a+)" \
                  << "\ncreates: curvature_local.dat (a+)" \
                  << "\ncreates: radius_curvature_global.dat (a+)" \
                  << "\ncreates: radius_curvature_local.dat (a+)" \
                  << std::endl;
        exit(1);
    }


    FILE * fp_proto_centroid_dist;
    fp_proto_centroid_dist = fopen("distance_proto_centroids.dat", "a+");

    // GLOBAL: global file ..
    FILE * fp_curvature_global;
    FILE * fp_radius_of_curvature_global;
    fp_curvature_global = fopen("curvature_global.dat", "a+");
    fp_radius_of_curvature_global = fopen("radius_curvature_global.dat", "a+");

    // LOCAL: local file ..
    FILE * fp_curvature_local;
    FILE * fp_radius_of_curvature_local;
    fp_curvature_local = fopen ("curvature_local.dat", "a+");
    fp_radius_of_curvature_local = fopen ("radius_curvature_local.dat", "a+");




    // data descriptions...
    // start frame
    fprintf(fp_proto_centroid_dist,"# Distance, tubulin monomer centroids, 44 Angstroms,typical.\n");
    fprintf(fp_proto_centroid_dist,"# start_frame: %d  stop: %d  step: %d\n",start,stop,step);

    fprintf(fp_curvature_global,"# Global Curvature: in Angstroms\n");
    fprintf(fp_curvature_global,"# start_frame: %d  stop: %d  step: %d\n",start,stop,step);

    fprintf(fp_radius_of_curvature_global,"# Global RadiusC:\n");
    fprintf(fp_radius_of_curvature_global,"# start_frame: %d  stop: %d  step: %d\n",start,stop,step);

    fprintf(fp_curvature_local,"# Local Curvature: (curv,left_z,right_z) for monomers 0-3,2-5,4-7,6-9,8-11.\n");
    fprintf(fp_curvature_local,"# start_frame: %d  stop: %d  step: %d\n",start,stop,step);

    fprintf(fp_radius_of_curvature_local,"# Local RadiusC: (RadiusC,left_z,right_z) for monomers 0-3,2-5,4-7,6-9,8-11.\n");
    fprintf(fp_radius_of_curvature_local,"# start_frame: %d  stop: %d  step: %d\n",start,stop,step);





    // << " <atom1> (i.e. 14)"        \
    // << " <atom2> (i.e. 257)" \

    // int resid1 = atoi(argv[6]);
    // int resid2 = atoi(argv[7]);



    /* ---------------------------------------------------------
       when I was using atoms closest to a point in space
       --------------------------------------------------------- */
    // // section 1
    // for ( int i=0; i <chains_to_use; i++) {
    //     chain_ref[i].ComputeCentroid();
    // }
    // Vector m;
    // int pos_low,pos_high;
    // for(int i=0; i<chains_to_use-1; i++){
    //     m = midpoint2(chain_ref[i].centroid,chain_ref[i+1].centroid);
    //     debug("reference_midpoint: %f %f %f\n",m.x,m.y,m.z);
    //     pos_high = get_pos_closest_to_point(chain_ref,i,m);
    //     pos_low = get_pos_closest_to_point(chain_ref,i+1,m);

    //     // switch to chain_later.
    //     chain_ref[i].pf_pos_high = pos_high;
    //     chain_later[i].pf_pos_high = pos_high;
    //     chain_ref[i+1].pf_pos_low = pos_low;
    //     chain_later[i+1].pf_pos_low = pos_low;
    // }

    /* ---------------------------------------------------------
       using command args to set atoms, 14,257
       --------------------------------------------------------- */
    // for(int i=0; i<chains_to_use; i++){
    //     chain_ref[i].pf_pos_high = resid2;
    //     chain_later[i].pf_pos_high = resid2;
    //     chain_ref[i+1].pf_pos_low = resid1;
    //     chain_later[i+1].pf_pos_low = resid1;
    // }


    // for ( int i=0;i<chains_to_use; i++) {
    //     printf("low_neighbor: %d (pos,index)  %d %d\n",i,chain_ref[i].pf_pos_low, \
    //            chain_ref[i].indices[chain_ref[i].pf_pos_low]);
    //     printf("high_neighbor: %d (pos,index)  %d %d\n",i,chain_ref[i].pf_pos_high, \
    //            chain_ref[i].indices[chain_ref[i].pf_pos_high]);
    // }



    Vector t0;
    Matrix Rz,Rx,Ry; // z x y


    // translate to origin
    t0 = translate_to_origin(chain_0,0,num_chains);
    // around Z to Y axis, less X.
    Rz = rotate_system_around_axis(chain_0,num_chains,2,1,\
                                   chain_0[0].centroid,\
                                   chain_0[chains_to_use-1].centroid,1);
    // around X to Y axis, less Z.
    Rx = rotate_system_around_axis(chain_0,num_chains,0,1,\
                                   chain_0[0].centroid,\
                                   chain_0[chains_to_use-1].centroid,1);
    // // print Centroids
    debug("<<<  ^^^after rotate_around_x_to_y  >>>\n");
    Vector proto_midpoint,proto_centroid;
    proto_centroid = get_centroid_from_centroids(chain_0,0,chains_to_use);
    proto_midpoint = midpoint2(chain_0[0].centroid,chain_0[chains_to_use-1].centroid);
    // around Y to
    Ry = rotate_system_around_axis(chain_0,num_chains,1,2,   \
                                   proto_midpoint,               \
                                   proto_centroid,1);


    // exit(0);
    // end of section 1.
#endif // CURVATURE_B

#ifdef ANGLE3CENTROID_B
    if (argc != 8) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference-PDB>" \
                  << " <Filename-timelater-DCD>" \
                  << " <num_chains>" \
                  << " <chains_ignore>" \
                  << " <start>"                  \
                  << " <stop>"                  \
                  << " <step>"                  \
                  << "\ncreates: mt_angles_ns.dat" \
                  << "\ncreates: mt_angles_ew.dat"
                  <<std::endl;
        // time ./run_segment_dcd_angle_mt mt.ref.pdb mt_d1_indent.dcd 209 1 0 500 25
        exit(1);
    }


    FILE * fp_angles_ns;
    FILE * fp_angles_ew;

    fp_angles_ns = fopen ("mt_angles_ns.dat", "a+");
    fp_angles_ew = fopen ("mt_angles_ew.dat", "a+");

    fprintf(fp_angles_ew, "# reference:\n");
    fprintf(fp_angles_ns, "# reference:\n");


    for( int i=0; i<chains_to_use; i++) {
        chain_ref[i].ComputeCentroid();
    }


    for( int i=0; i<chains_to_use; i++) {
        get_latangle_from35centroids(chain_ref,chain_ref,i,fp_angles_ns,0); // N,S
        get_latangle_from35centroids(chain_ref,chain_ref,i,fp_angles_ew,1); // E,W
    }

    for( int i=0; i<chains_to_use; i++) {
        chain_later[i].ComputeCentroid();
    }

    // printf("\n\ngetting latangle from chain_ref (REFERENCE)!\n");
    // printf("entering ANGLE3CENTROID!!!!\n");
    // debug("\n\ngetting latangle from chain_later!\n\n");

    // start or frame_position
    fprintf(fp_angles_ns, "\n# time later:   C|South|North (frame_start:%d)\n",start);
    fprintf(fp_angles_ns, "# if proto: 102(3),324(11),546(19),768(27),9810(35).");

    fprintf(fp_angles_ew, "\n# time later:   C|West|East   (frame_start:%d)\n",start);
#endif // ANGLE3CENTROID_B




#if defined(CONTACTMAP_SM1) || defined(CONTACTMAP_MT1)
    // map
    if (argc != 8) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference-PDB>" \
                  << " <Filename-timelater-DCD>" \
                  << " <num_chains>" \
                  << " <chains_ignore>" \
                  << " <start>" \
                  << " <stop>" \
                  << " <step>" \
                  << "\nfiles created: contacts_{n,e,s,w,intra(1)}.dat" \
                  <<std::endl;
        // << "\ncase: (intra) if num_chains = 1, chains_ignore = 0" \
        exit(1);
    }

    std::cout << "\n\nWelcome to contact map analysis!\n\n\n" \
              << "the contact definition involves a 8A-C_ALPHA atom spherical cutoff\n"
              << "the persistence criterion involves a <8A or less than original distance + 2.0 A\n"
              << "files to be created include contacts_{n,e,s,w,intra(1)}.dat.\n\n" \
              <<std::endl;

#endif // #if defined(CONTACTMAP_SM1) || defined(CONTACTMAP_MT1)

#ifdef CONTACTMAP_SM1
    // for (int i=0; i<chains_to_use; i++) {
    // for (int i=0; i<3; i++) {
    Contact map_ref[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];

    topo_count_map(map_ref);
    topo_build_intra(chain_ref,0,map_ref,8.0); // *chain,int chain_num,Contact (*map)[MAX_CONTACTS],cutoff
    memcpy(chain_ref[0].contacts,map_ref, sizeof(chain_ref[0].contacts));
    topo_count_map(map_ref);
    topo_clear_map(map_ref);
    // }

    FILE * fp_contacts;
    fp_contacts = fopen("contacts_intra.dat", "w+");

#ifdef CONTACTSBYRESFILE1
    FILE * fp_contactsresfile1;
    fp_contactsresfile1 = fopen("contacts_by_residue.dat","w+");

#endif // CONTACTSBYRESFILE1

#endif // CONTACTMAP_SM1

#ifdef CONTACTMAP_MT1
    Contact map_ref[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];

    int num_contacts;
    num_contacts = 0;

    // printf("done with intra.\n");
    // for (int i=0; i<3; i++) {
    //     count_map(chain_ref[i].contacts);
    // }

    // intercontacts, with 8.0 A cutoff
    // Contact map_ref[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];


    // ONLY going to run if System has contactsLonN, LatE, LonS, and LatW defined.
    // if(chains_to_use > 1){
    for (int i=0; i<chains_to_use; i++) {
        // for (int i=0; i<1; i++) {


        // N
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LonN,map_ref,8.0);
        memcpy(chain_ref[i].contactsLonN,map_ref, sizeof(chain_ref[i].contactsLonN));
        // topo_count_map(map_ref);
        num_contacts = topo_count_map(map_ref);
#ifdef DIMERMAP
        maprow_contacts.push_back(num_contacts);
#endif
        topo_clear_map(map_ref);

        // E
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LatE,map_ref,8.0);
        memcpy(chain_ref[i].contactsLatE,map_ref, sizeof(chain_ref[i].contactsLatE));
        // topo_count_map(map_ref);
        num_contacts = topo_count_map(map_ref);
#ifdef DIMERMAP
        maprow_contacts.push_back(num_contacts);
#endif
        topo_clear_map(map_ref);

        // S
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LonS,map_ref,8.0);
        memcpy(chain_ref[i].contactsLonS,map_ref, sizeof(chain_ref[i].contactsLonS));
        // topo_count_map(map_ref);
        num_contacts = topo_count_map(map_ref);
#ifdef DIMERMAP
        maprow_contacts.push_back(num_contacts);
#endif
        topo_clear_map(map_ref);

        // W
        topo_count_map(map_ref);
        topo_build_inter(chain_ref,i,chain_ref[i].LatW,map_ref,8.0);
        memcpy(chain_ref[i].contactsLatW,map_ref, sizeof(chain_ref[i].contactsLatW));
        // topo_count_map(map_ref);
        num_contacts = topo_count_map(map_ref);
#ifdef DIMERMAP
        maprow_contacts.push_back(num_contacts);
#endif
        topo_clear_map(map_ref);


#ifdef DIMERMAP
        map_contacts.insert( std::pair<int, std::vector<int> >(i,maprow_contacts));
        maprow_contacts.clear();
#endif
    }



#ifdef CONFILE
    FILE * fp_contacts_n;
    fp_contacts_n = fopen("contacts_n.dat", "w+");
    FILE * fp_contacts_e;
    fp_contacts_e = fopen("contacts_e.dat", "w+");
    FILE * fp_contacts_s;
    fp_contacts_s = fopen("contacts_s.dat", "w+");
    FILE * fp_contacts_w;
    fp_contacts_w = fopen("contacts_w.dat", "w+");
#endif // CONFILE


    // OPEN FILE: assemble dimer contact counts.
    FILE * fp_dimercontacts;
    fp_dimercontacts = fopen("contacts_by_dimer.dat","w+");
    fprintf(fp_dimercontacts,"# start_frame: %d  stop: %d  step: %d\n",start,stop,step);
    fprintf(fp_dimercontacts,"# 0-5  alpha index findex | beta index findex\n");
    fprintf(fp_dimercontacts,"# 6-10 SINEW: southern, internal, northern, eastern, western\n");

    std::vector<int> dimers_alpha;
    for (int i=0; i<chains_to_use; i++){

        // alpha-tubulin
        if (chain_ref[i].num_atoms_ca >= 430){
            dimers_alpha.push_back(i);
        };

    }



    // assemble dimer contact counts.
    for (int j=0; j<dimers_alpha.size(); j++) {
        debug("%d ---> %d\n",j,dimers_alpha[j]);
        topo_sort_dimer_contacts_initially(chain_ref,dimers_alpha[j],fp_dimercontacts);
    }


    // Map cm1;
    // for (int i=0; i<chains_to_use; i++){

    //     cm1.build_map_default(chain_ref,i);
    // }
    // exit(0);

#endif // CONTACTMAP_MT1


    /* ---------------------------------------------------------
       Microtubules (1 of 2)
       --------------------------------------------------------- */
#ifdef MTCON1
    std::cout << "this segment of code is fully deprecated." << chains_to_use << std::endl;
    exit(1);


    std::cout << "MT! using chains: " << chains_to_use << std::endl;


    Vector centroid_front;
    Vector centroid_back;
    Vector mt_axis;
    Vector mt_axis_norm;
    // centroid_front.x = 0.0, centroid_front.y = 0.0, centroid_front.z = 0.0;
    // std::cout << centroid_front.x << centroid_front.y << centroid_front.z << std::endl;
    // std::cout << centroid_back.x << centroid_back.y << centroid_back.z << std::endl;
    // std::cout << mt_axis.x << mt_axis.y << mt_axis.z << std::endl;
    // std::cout << mt_axis_norm.x << mt_axis_norm.y << mt_axis_norm.z << std::endl;


    // ComputeCentroid()
    for ( int i=0; i < num_chains; i++) {
        chain_ref[i].ComputeCentroid();
        // EXACTLY THE SAME!
        // printf("%f %f %f\n",chain_ref[i].centroid.x, \
        //        chain_ref[i].centroid.y,\
        //        chain_ref[i].centroid.z);

        // printf("%f %f %f\n",&ref_centroid_pos[i].x,\
        //        &ref_centroid_pos[i].y,\
        //        &ref_centroid_pos[i].z);
    }


    // // Array of Centroids (1 way)
    // Vector ref_centroid_pos[chains_to_use-1];
    // // std::vector<Vector> ref_centroid_pos[chains_to_use-1];
    // for ( int i=0; i<chains_to_use; i ++ ) {
    //     get_tubulin_centroid(&chain_ref[i],&ref_centroid_pos[i]);
    //     // get_tubulin_centroid(&chain_ref[i],&ref_centroid_pos[i]);
    // }
    // // get_centroid_from_array(ref_centroid_pos, 0, chains_to_use/2, &centroid_front);
    // // get_centroid_from_array(ref_centroid_pos, chains_to_use/2 + 1, chains_to_use, &centroid_back);
    // get_centroid_from_array(ref_centroid_pos, 0, 30, &centroid_front);
    // get_centroid_from_array(ref_centroid_pos, chains_to_use - 30, chains_to_use, &centroid_back);

    // // print centroids (1)
    // get_vector(centroid_front, centroid_back, &mt_axis);
    // normalize(mt_axis, &mt_axis_norm);
    // debug("Centroid_front: %f %f %f",centroid_front.x,centroid_front.y,centroid_front.z);
    // debug("Centroid_back: %f %f %f",centroid_back.x,centroid_back.y,centroid_back.z);
    // debug("mt_axis: %f %f %f",mt_axis.x,mt_axis.y,mt_axis.z);
    // debug("mt_axis_norm: %f %f %f",mt_axis_norm.x,mt_axis_norm.y,mt_axis_norm.z);



    // centroid from Centroids (2nd way)
    // centroid_front = get_centroid_from_centroids(chain_ref,0, chains_to_use/2);
    // centroid_back = get_centroid_from_centroids(chain_ref,chains_to_use/2 + 1, chains_to_use);
    centroid_front = get_centroid_from_centroids(chain_ref, 0, 40);
    centroid_back = get_centroid_from_centroids(chain_ref, chains_to_use - 46, chains_to_use-15);


    // print centroids (2)
    get_vector(centroid_front, centroid_back, &mt_axis);
    normalize(mt_axis, &mt_axis_norm);
    debug("Centroid_front: %f %f %f\n",centroid_front.x,centroid_front.y,centroid_front.z);
    debug("Centroid_back: %f %f %f\n",centroid_back.x,centroid_back.y,centroid_back.z);
    debug("mt_axis: %f %f %f\n",mt_axis.x,mt_axis.y,mt_axis.z);
    debug("mt_axis_norm: %f %f %f\n",mt_axis_norm.x,mt_axis_norm.y,mt_axis_norm.z);


    printf(".. extending analysis to contact interfaces .......\n");
    compute_nearest_neighbors(chain_ref,chains_to_use);

    printf(".. determining latitudinal and longitudinal neighbors .......\n");
    determine_latlon_neighbors(chain_ref,chains_to_use,mt_axis_norm);

    for ( int i=0; i < num_chains; i++) {
    // for ( int i=5; i < 15; i++) {
        debug("chain:(%d) \n",i);
        debug("\n%3d >> N:%3d S:%3d W:%3d E:%3d\n",i,chain_ref[i].LonN,chain_ref[i].LonS,chain_ref[i].LatW,chain_ref[i].LatE);
    }

    for ( int i=0; i < chains_to_use; i++) {
        // printf("chain:(%d) ",i);
        // printf("%3d >> N:%3d S:%3d W:%3d E:%3d\n",i,chain_ref[i].LonN,chain_ref[i].LonS,chain_ref[i].LatW,chain_ref[i].LatE);
        // if ( i%20 == 0 ) {
        //     printf("%d   %f %f %f\n",i,chain_ref[i].centroid.x,chain_ref[i].centroid.y,chain_ref[i].centroid.z);
        //     printf("%d   %f %f %f\n",i,ref_centroid_pos[i].x,ref_centroid_pos[i].y,ref_centroid_pos[i].z);
        // }
        if ( chain_ref[i].LonN == -1 ) {
            printf("no N  %d:%d [%d|%d] %d\n",i,chain_ref[i].LonN,chain_ref[i].index,chain_ref[i].findex,\
                   chain_ref[i].num_atoms_ca);
        }
        if ( chain_ref[i].LatE == -1 ) {
            printf("no E  %d:%d [%d|%d] %d\n",i,chain_ref[i].LatE,chain_ref[i].index,chain_ref[i].findex,\
                   chain_ref[i].num_atoms_ca);
        }
        if ( chain_ref[i].LonS == -1 ) {
            printf("no S  %d:%d [%d|%d] %d\n",i,chain_ref[i].LonS,chain_ref[i].index,chain_ref[i].findex,\
                   chain_ref[i].num_atoms_ca);
        }
        if ( chain_ref[i].LatW == -1 ) {
            printf("no W  %d:%d [%d|%d] %d\n",i,chain_ref[i].LatW,chain_ref[i].index,chain_ref[i].findex,\
                   chain_ref[i].num_atoms_ca);
        }
    }
    // exit(0);


    printf(".. getting latitudinal and longitudinal contacts .......\n");
    // N E S W | 1 2 3 4

    for ( int i=0; i<chains_to_use; i++ ) {
    // for ( int i=33; i<34; i++ ) {

// #undef MOLECULE_INITIAL_SIZE
// #define MOLECULE_INITIAL_SIZE chain_ref[i].num_atoms_ca


        // PRINT HERE
        // print_contact_map(chain_ref[i].contactsLonN); // 0
        // printf("before_counter: %d\n",chain_ref[i].LonN_total2x);
        Contact map1[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter_monomer(chain_ref,i,map1,1);
        memcpy(chain_ref[i].contactsLonN, map1, sizeof(chain_ref[i].contactsLonN));
        // printf("after_counter:(N) %d\n",chain_ref[i].LonN_total2x);
        // print_contact_map(chain_ref[i].contactsLonN); // 0
        // exit(0);

        Contact map2[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter_monomer(chain_ref,i,map2,2);
        memcpy(chain_ref[i].contactsLatE, map2, sizeof(chain_ref[i].contactsLatE));
        // printf("after_counter:(E) %d\n",chain_ref[i].LatE_total2x);
        Contact map3[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter_monomer(chain_ref,i,map3,3);
        memcpy(chain_ref[i].contactsLonS, map3, sizeof(chain_ref[i].contactsLonS));
        // printf("after_counter:(S) %d\n",chain_ref[i].LonS_total2x);
        Contact map4[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter_monomer(chain_ref,i,map4,4);
        memcpy(chain_ref[i].contactsLatW, map4, sizeof(chain_ref[i].contactsLatW));
        // printf("after_counter:(W) %d\n",chain_ref[i].LatW_total2x);
        // exit(0);



        // // PRINT HERE
        // // print_contact_map(chain_ref[i].contactsLonN); // 0
        // // printf("before_counter: %d\n",chain_ref[i].LonN_total2x);
        // Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
        // build_contact_map_inter_monomer(chain_ref,i,map,1);
        // memcpy(chain_ref[i].contactsLonN, map, sizeof(chain_ref[i].contactsLonN));
        // // printf("after_counter: %d\n",chain_ref[i].LonN_total2x);
        // // print_contact_map(chain_ref[i].contactsLonN); // 0
        // // exit(0);


        // // printf("before_counter: %d\n",chain_ref[i].LatE_total2x);
        // Contact map2[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
        // // build_contact_map_inter_monomer(chain_ref,i,map2,2);
        // build_contact_map_inter_monomer(chain_ref,i,map2,2);
        // memcpy(chain_ref[i].contactsLatE, map2, sizeof(chain_ref[i].contactsLatE));
        // // printf("after_counter: %d\n",chain_ref[i].LatE_total2x);
        // // print_contact_map(chain_ref[i].contactsLatE); // 0

        // // printf("before_counter: %d\n",chain_ref[i].LonS_total2x);
        // Contact map3[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
        // build_contact_map_inter_monomer(chain_ref,i,map3,3);
        // memcpy(chain_ref[i].contactsLonS, map3, sizeof(chain_ref[i].contactsLonS));
        // // printf("after_counter: %d\n",chain_ref[i].LonS_total2x);
        // // print_contact_map(chain_ref[i].contactsLonS); // 0

        // // printf("before_counter: %d\n",chain_ref[i].LatW_total2x);
        // Contact map4[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
        // build_contact_map_inter_monomer(chain_ref,i,map4,4);
        // memcpy(chain_ref[i].contactsLatW, map4, sizeof(chain_ref[i].contactsLatW));
        // // printf("after_counter: %d\n",chain_ref[i].LatW_total2x);
        // // print_contact_map(chain_ref[i].contactsLatW); // 0
    }
    // exit(0);


// #ifdef MTTOP
//     // for ( int i=0; i<chains_to_use; i++ ) {
//     //     // PRINT HERE
//     //     // print_contact_map(chain_ref[i].contactsLonN); // 0
//     //     // printf("before_counter: %d\n",chain_ref[i].LonN_total2x);
//     //     Contact map1[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
//     //     build_contact_map_inter_monomer(chain_ref,i,map1,1);
//     //     memcpy(chain_ref[i].contactsLonN, map1, sizeof(chain_ref[i].contactsLonN));
//     // }

//     printf("building intra __contacts__ map..\n");
//     for ( int i=0; i<chains_to_use; i++ ) {
//         debug("chain:(%d) ",i);
//         build_contact_map(&chain_ref[i]);
//     }

//     printf("assigning eh..\n");
//     for ( int i=0; i<chains_to_use; i++ ) {
//         debug("chain:(%d) ",i);
//         assign_contact_eh(chain_ref[i].contacts,1.33333); // gets the rest of the 0.0;
//         assign_contact_eh(chain_ref[i].contactsLonN,1.99999);
//         assign_contact_eh(chain_ref[i].contactsLatW,1.88888);
//     }

//     printf("printing topology..\n");
//     for ( int i=0; i<chains_to_use; i++ ) {
//         fprintf_contact_map_gsop(chain_ref[i].contacts);
//         fprintf_contact_map_gsop(chain_ref[i].contactsLonN);
//         fprintf_contact_map_gsop(chain_ref[i].contactsLatW);
//     }


// #endif // MTTOP

#endif // MTCON1 Section 1. END



#ifdef DIMERMAP_1B
    printf("starting section dimermap_1B now ..\n");

#ifndef NDEBUG
    // print map, check.
    printf("contacts:\n");
    for (int j=0; j<map_contacts.size(); j++)
    {
        printf("\t%d: %d",j,map_contacts[j].size());
        for (int k=0; k<map_contacts[j].size(); k++)
        {
            printf("  %d-%d",k,map_contacts[j][k]);
        }
        printf("\n");
    }
#endif // NDEBUG


    // angles:
    for (int i=0; i<dimermap_alphas.size(); i++)
    {
        //                                                 *2nd is "chain_later"
        maprow_angles = get_3angles_dimerbyalpha(chain_ref,chain_ref,dimermap_alphas[i]);
        map_angles.insert( std::pair<int, std::vector<double> >(i,maprow_angles));
        maprow_angles.clear();
    }


    // curvature:
    int a_i, b_i;
    Curvature curv_a, curv_b, curv_lon;

    for (int i=0; i<chains_to_use; i++ )
    {
        chain_ref[i].ComputeCentroid();
    }

    for (int i=0; i<dimermap_alphas.size(); i++)
    {
        // alat
        a_i = dimermap_alphas[i];
        curv_a.p1 = chain_ref[a_i].centroid;
        curv_a.p2 = chain_ref[chain_ref[a_i].LatW].centroid;
        curv_a.p3 = chain_ref[chain_ref[a_i].LatE].centroid;
        // curv_a.print_points();
        // printf("\t-->\n");
        // curv_a.points_ascending_z();
        // curv_a.print_points();
        // curv_a.get_circumscribed_circle();
        curv_a.get_radius_circumscribed();
        debug(">: %f %f\n",curv_a.radius,curv_a.curvature);
        // printf("alat ^\n");
        maprow_curv.push_back(curv_a.radius);

        // blat
        b_i = chain_ref[a_i].LonN;
        curv_b.p1 = chain_ref[b_i].centroid;
        curv_b.p2 = chain_ref[chain_ref[b_i].LatW].centroid;
        curv_b.p3 = chain_ref[chain_ref[b_i].LatE].centroid;
        // curv_b.print_points();
        // printf("\t-->\n");
        // curv_b.points_ascending_z();
        // curv_b.print_points();
        // curv_b.get_circumscribed_circle();
        curv_b.get_radius_circumscribed();
        // printf(">: %f %f\n",curv_b.radius,curv_b.curvature);
        // printf("blat ^\n");
        maprow_curv.push_back(curv_b.radius);

        // blon
        curv_lon.p1 = chain_ref[a_i].centroid;
        curv_lon.p2 = chain_ref[b_i].centroid;
        curv_lon.p3 = chain_ref[chain_ref[b_i].LonN].centroid;
        // curv_lon.print_points();
        // printf("\t-->\n");
        // curv_lon.points_ascending_z();
        // curv_lon.print_points();
        // curv_lon.get_circumscribed_circle();
        curv_lon.get_radius_circumscribed();
        // printf(">: %f %f\n",curv_lon.radius,curv_lon.curvature);
        // printf("blon ^\n");
        maprow_curv.push_back(curv_lon.radius);

        map_curv.insert(std::pair<int, std::vector<double> >(i,maprow_curv));
        // map_contacts.insert( std::pair<int, std::vector<int> >(i,maprow_contacts));

        maprow_curv.clear();
    }
    // exit(0

#ifndef NDEBUG
    // print map, check.
    printf("angles:\n");
    for (int j=0; j<map_angles.size(); j++)
    {
        printf("\t%d: %d",j,map_angles[j].size());
        for (int k=0; k<map_angles[j].size(); k++)
        {
            printf("   %d %f",k,map_angles[j][k]);
        }
        printf("\n");
    }
#endif // NDEBUG


    // max:
    fprintf(fp_mt_analysis,"# max/begin:\n");

    // write data!
    // counts of each contacts
    int dimer_east,dimer_west,dimer_ext,alpha,beta;
    dimer_east = dimer_west = dimer_ext = alpha = beta = 0;

    for (int i=0; i<dimermap_alphas.size(); i++)
    {
        // sum contacts
        alpha = dimermap_alphas[i];
        beta = chain_ref[alpha].LonN;
        dimer_east = map_contacts[alpha][1] + map_contacts[beta][1];
        dimer_west = map_contacts[alpha][3] + map_contacts[beta][3];
        dimer_ext = dimer_east + dimer_west + map_contacts[alpha][2] + map_contacts[beta][0];


        // WRITE DATA!
        // S,I,N,E,W,E
        fprintf(fp_mt_analysis,"%3d %6d %6d %3d %6d %6d  %3d %3d %3d %3d %3d %3d", \
                alpha,chain_ref[alpha].index,chain_ref[alpha].findex,   \
                beta,chain_ref[beta].index,chain_ref[beta].findex,      \
                map_contacts[alpha][2],map_contacts[alpha][0],          \
                map_contacts[beta][0],dimer_east,dimer_west,dimer_ext);

        // alpha lateral, beta lateral, beta longitudinal
        fprintf(fp_mt_analysis," %5.1f %5.1f %5.1f",  \
                map_angles[i][0],map_angles[i][1],map_angles[i][2]);

        // alat, blat, blon
        fprintf(fp_mt_analysis,"  %7.1f %7.1f %7.1f\n", \
                map_curv[i][0], \
                map_curv[i][1], \
                map_curv[i][2]);


    }
    // fprintf(fp_mt_analysis,"\n");


    // map clear.
    map_contacts.clear();
    map_angles.clear();
    map_curv.clear();

#endif // DIMERMAP_1B





    /* ---------------------------------------------------------
       DCD Preface.
       --------------------------------------------------------- */



    // Start DCD
    molfile_timestep_t timestep;
    void *v;
    dcdhandle *dcd;
    // int i, natoms;
    int natoms;
    float sizeMB =0.0, totalMB = 0.0;
    double starttime, endtime, totaltime = 0.0;

    // // 1
    // while (--argc) {
    //     ++argv;
    //     natoms = 0;
    //     v = open_dcd_read(*argv, "dcd", &natoms);
    //     if (!v) {
    //         fprintf(stderr, "main) open_dcd_read failed for file %s\n", *argv);
    //         return 1;
    //     }
    //     dcd = (dcdhandle *)v;
    //     sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    //     totalMB += sizeMB;
    //     printf("main) file: %s\n", *argv);
    //     printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);

    //     // starttime = time_of_day();
    //     timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    //     for (i=0; i<dcd->nsets; i++) {
    //         int rc = read_next_timestep(v, natoms, &timestep);
    //         if (rc) {
    //             fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
    //             return 1;
    //         }
    //     }
    //     // endtime = time_of_day();
    //     close_file_read(v);
    //     // totaltime += endtime - starttime;
    //     // printf("  Time: %5.1f seconds\n", endtime - starttime);
    //     // printf("  Speed: %5.1f MB/sec, %5.1f timesteps/sec\n", sizeMB \
    //     //        / (endtime - starttime), (dcd->nsets / (endtime - starttime)));
    // }
    // printf("Overall Size: %6.1f MB\n", totalMB);
    // // printf("Overall Time: %6.1f seconds\n", totaltime);
    // // printf("Overall Speed: %5.1f MB/sec\n", totalMB / totaltime);



    // int atoms_in_chain;
    // 2. to read a dcd.
    natoms = 0;
    v = open_dcd_read(argv[2], "dcd", &natoms);
    if (!v) {
        fprintf(stderr, "main) open_dcd_read failed for file %s\n", *argv);
        return 1;
    }
    dcd = (dcdhandle *)v;
    sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
    totalMB += sizeMB;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    // timestep.velocities = (float *)malloc(3*sizeof(float)*natoms);

    // printf("main) file: %s\n", *argv);
    // printf("  %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);


    // close_file_read(v);
    // END DCD
    // typedef struct {
    //     fio_fd fd;
    //     int natoms; 382
    //     int nsets; 17501 (0-17500)
    //     int setsread;
    //     int istart;
    //     int nsavc;
    //     double delta;
    //     int nfixed;
    //     float *x, *y, *z; ->x[0-381];
    //     int *freeind;
    //     float *fixedcoords;
    //     int reverse;
    //     int charmm;
    //     int first;
    //     int with_unitcell;
    // } dcdhandle;
    printf("--------START HERE-------------\n");

    // dcd
    // 0: (pdb) | ref | chain_0 (from dcd) | chain_later
    // int advance_dcd(int numframes,int frame,dcdhandle *v,int natoms,molfile_timestep_t *timestep);

    // printf("ref: %f\n",chain_ref[0].pos[105].y);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);

    // frame_position = advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 1
    // load_dcd_to_chain(dcd,chain_0,num_chains);
    // printf("frame_position(0): %d\n",frame_position);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);

    // frame_position = advance_dcd(dcd->nsets,1,dcd,natoms,&timestep); // 1
    // load_dcd_to_chain(dcd,chain_0,num_chains);
    // printf("frame_position: %d\n",frame_position);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);

    // frame_position = advance_dcd(dcd->nsets,10,dcd,natoms,&timestep); // 1
    // load_dcd_to_chain(dcd,chain_0,num_chains);
    // printf("frame_position: %d\n",frame_position);
    // printf("0-105: %f\n",chain_0[0].pos[105].y);


    // printf("beginning new loop.\n");
    // int countd;

    // countd = 0;
    // for(int df=start; df<=stop; df+=step){
    //     countd +=1;
    //     // dummy = my_dcd_read(df);
    //     frame_position = advance_dcd(dcd->nsets,df,dcd,natoms,&timestep);
    //     printf("frame_position:--->  %d  <--- %d %d\n",frame_position,countd,df);
    //     load_dcd_to_chain(dcd,chain_later,num_chains);
    //     printf("0-105: %f\n",chain_later[0].pos[105].y);

    // }
    // exit(0);

    // nil.
    printf("ref-findex(%d): %f\n",chain_ref[0].findex,chain_ref[0].pos[chain_ref[0].findex].y);
    printf("0-findex(%d): %f\n",chain_0[0].findex,chain_0[0].pos[chain_0[0].findex].y);
    printf("later-findex(%d): %f\n",chain_later[0].findex,chain_later[0].pos[chain_later[0].findex].y);

    frame_position = 1;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 1st advance. 1-vmd
    load_dcd_to_chain(dcd,chain_0,num_chains);
    printf("frame_position: %d\n",frame_position);

    printf("ref-findex(%d): %f\n",chain_ref[0].findex,chain_ref[0].pos[chain_ref[0].findex].y);
    printf("0-findex(%d): %f\n",chain_0[0].findex,chain_0[0].pos[chain_0[0].findex].y);
    printf("later-findex(%d): %f\n",chain_later[0].findex,chain_later[0].pos[chain_later[0].findex].y);

    double pos1,pos2;
    pos1 = chain_0[0].pos[chain_0[0].findex].y;
    pos2 = chain_later[0].pos[chain_later[0].findex].y;
    if(pos2 - pos1 < 0.0001){
        printf("Your reference and chain_0 are likely the same.\n");
        printf("You are likely using the PDB_COORDS from which the simulation began.\n");
        sleep(1);
    } else {
        printf("WARNING: your reference and chain_0 are different.\n");
        printf("You may be using PDB_COORDS from which the simulation was not begun.\n");
        sleep(2);
    }

    frame_position = 2;
    advance_dcd(dcd->nsets,0,dcd,natoms,&timestep); // 2nd. 2-vmd
    load_dcd_to_chain(dcd,chain_later,num_chains);
    printf("frame_position: %d\n",frame_position);
    // exit(0);

    // Advancing Rules.
    // ----------------
    // example. step size -> 5.
    // int advance_size = atoi(argv[3]) - 1;

    // step.
    int advance_size = step - 1; // 0 counts, so advance_size of 4, advances by 5.

    // stop.
    if(stop > dcd->nsets){
        stop = dcd->nsets;
        printf("use stop value: %d\n",stop);
    }

    // // start.
    // debug("starting frame: %d\n",start);
    // if((start > 2) && (step < start)) {
    //     // for (int nset1=2; nset1<start; nset1 += advance_size + 1 ) {
    //     for (int nset1=2; nset1<start-step; nset1 += 1 ) {
    //         // for (int nset1=2; nset1<dcd->nsets; nset1 += step + 1) {
    //         debug("forwarding --> current: %d\n",nset1);
    //         frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);
    //         debug("forwarding --> frame_position: %d\n",frame_position);

    //         // if ( nset1 + advance_size >= start) {
    //         //     // if ( nset1 + step >= dcd->nsets ) {
    //         //     break;
    //         // }
    //         // frame_position is returned value;
    //         // frame_position += advance_dcd(dcd->nsets,step,dcd,natoms,&timestep);
    //         // load_dcd_to_chain(dcd,chain_later,num_chains);
    //     }
    // } else if ((start > 2) && (step > start - 2)) {
    //     for (int nset1=2; nset1<start; nset1 += 1 ) {
    //         debug("forwarding --> current: %d\n",nset1);
    //         frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);
    //         debug("forwarding --> frame_position: %d\n",frame_position);
    //     }
    // }

    // sleep(0.5);
    // exit(0);


    for (int nset1=2; nset1<start; nset1 += 1 ) {
        // for (int nset1=2; nset1<dcd->nsets; nset1 += step + 1) {
        // debug("forwarding --> current: %d\n",nset1);
        frame_position += advance_dcd(dcd->nsets,0,dcd,natoms,&timestep);

        load_dcd_to_chain(dcd,chain_later,num_chains);

        debug("forwarding --> frame_position: %d\n",frame_position);
    }

    // Get initial starting point.
    printf("--> fast forwarded. to frame: %d\n",frame_position);


    /* ---------------------------------------------------------
       major for loop begin || doloop.
       --------------------------------------------------------- */
    // for (int nset2=frame_position; nset2<stop; nset2 += advance_size + 1) {
    //     // for (int nset=frame_position; nset<dcd->nsets; nset += advance_size + 1 ) {

    //     if ( nset2 + advance_size >= stop) {
    //         // if ( nset + advance_size >= dcd->nsets ) {
    //         break;
    //     }

    //     frame_position += advance_dcd(dcd->nsets,advance_size,dcd,natoms,&timestep);
    //     load_dcd_to_chain(dcd,chain_later,num_chains);
    //     debug("current: %d\n",nset2);
    //     printf("--> frame_position: %d\n",frame_position);
    //     // check.
    //     // printf("ref-findex(%d): %f\n",chain_ref[0].findex,chain_ref[0].pos[chain_ref[0].findex].y);
    //     // printf("0-findex(%d): %f\n",chain_0[0].findex,chain_0[0].pos[chain_0[0].findex].y);
    //     // printf("later-findex(%d): %f\n",chain_later[0].findex,chain_later[0].pos[chain_later[0].findex].y);
    //     // continue;

    int nset2;
    nset2 = frame_position;
    do {

#endif //DCDREAD



#ifdef DEFAULT
        if (argc != 5) {
            std::cout << "Usage: " << argv[0] \
                      << " <Filename-reference>" \
                      << " <Filename-timelater>" \
                      << " ... ?opt args .." \
                      <<std::endl;
            exit(1);
        }
    // run_segment  ref   timelater maxchains ignorechains
    // argc 1        2          3        4           5
    // argv 0        1          2        3           4
    num_chains = atoi(argv[3]); // argv = 3, chains are 0,1,2
    chains_ignore = atoi(argv[4]);
    chains_to_use = num_chains - chains_ignore;



    // REFERENCE
    Chain *chain_ref;
    try {
        chain_ref = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }
    // int i=0;
    // int num_chain_ref=0, chains_to_useref=0, chains_ignoreref=0;


    for ( int i = 0; i < num_chains; i ++ ) {
        debug("building chain --> %d\n",i);
        strcpy(chain_ref[i].filename,argv[1]);
        chain_ref[i].chainid = i;
        if ( i == 0 ) {
            debug("building first chain (%d)\n", i);
            chain_ref[i].index = 0;
            chain_ref[i].file_line_begin = 0;
        } else {

            debug(" building later chain --> %d\n",i);
            chain_ref[i].index = chain_ref[i-1].findex + 1;
            chain_ref[i].file_line_begin = chain_ref[i-1].file_line_end; // or +1?
        }
        // Read PDB File.
        ReadLines(&chain_ref[i]);
        chain_ref[i].findex = chain_ref[i].index + chain_ref[i].num_atoms_ca - 1;
        chain_ref[i].assign_indices();


#ifdef INFO
        if (chain_ref[i].num_atoms_ca > 0) {chain_ref[i].print_prop();}
#endif // INFO
    }
    std::cout << " <<<<<<<< REFERENCE has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;


    Chain *chain_later;
    try {
        chain_later = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }

    for ( int i = 0; i < num_chains; i ++ ) {
        debug("building chain --> %d\n",i);
        strcpy(chain_later[i].filename,argv[2]);
        chain_later[i].chainid = i;
        if ( i == 0 ) {

            debug("building first chain (%d)\n", i);
            chain_later[i].index = 0;
            chain_later[i].file_line_begin = 0;
        } else {

            debug(" building later chain --> %d\n",i);
            chain_later[i].index = chain_later[i-1].findex + 1;
            chain_later[i].file_line_begin = chain_later[i-1].file_line_end; // or +1?
        }
        // Read PDB File.
        ReadLines(&chain_later[i]);
        chain_later[i].findex = chain_later[i].index + chain_later[i].num_atoms_ca - 1;
        chain_later[i].assign_indices();

#ifdef INFO
        if (chain_ref[i].num_atoms_ca > 0) {chain_ref[i].print_prop();}
#endif //
        // chain_later[i].print_prop();
    }
    std::cout << " <<<<<<<< TIME_LATER_STRUCT has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;


    // Acquiring Done.
    std::cout << "all chains acquired.\n" << std::endl;
#endif // DEFAULT


#ifdef ONEMOL
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference>" \
                  << " <num_chains>" \
                  << " <chain_ignore>" \
                  <<std::endl;
        exit(1);
    }
    // run_segment timelater maxchains ignorechains
    // argc 1        2          3        4
    // argv 0        1          2        3           4
    num_chains = atoi(argv[2]); // argv = 3, chains are 0,1,2
    chains_ignore = atoi(argv[3]);
    chains_to_use = num_chains - chains_ignore;


    Chain *chain_ref;
    try {
        chain_ref = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }
    // int i=0;
    // int num_chain_ref=0, chains_to_useref=0, chains_ignoreref=0;
    // int num_chain_ref=0, chains_to_useref=0, chains_ignoreref=0;

    for ( int i = 0; i < num_chains; i ++ ) {
        debug("building chain --> %d\n",i);
        strcpy(chain_ref[i].filename,argv[1]);
        chain_ref[i].chainid = i;
        if ( i == 0 ) {
            debug("building first chain (%d)\n", i);
            chain_ref[i].index = 0;
            chain_ref[i].file_line_begin = 0;
        } else {
            debug(" building later chain --> %d\n",i);
            chain_ref[i].index = chain_ref[i-1].findex + 1;
            chain_ref[i].file_line_begin = chain_ref[i-1].file_line_end; // or +1?
        }

        // Read PDB File.
        ReadLines(&chain_ref[i]);
        chain_ref[i].findex = chain_ref[i].index + chain_ref[i].num_atoms_ca - 1;
        chain_ref[i].assign_indices();
        // chain_ref[i].print_prop();

#ifdef INFO
        if (chain_ref[i].num_atoms_ca > 0) {chain_ref[i].print_prop();}
#endif // INFO
    }
    std::cout << " <<<<<<<< ONE structure has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;
#endif // ONEMOL

#ifdef ALLATOM
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-reference>" \
                  << " <num_chains>" \
                  << " <chain_ignore>" \
                  <<std::endl;
        exit(1);
    }
    // ./run_segment PDB maxchains ignorechains
    // argc 1        2          3        4
    // argv 0        1          2        3           4
    num_chains = atoi(argv[2]); // argv = 3, chains are 0,1,2
    chains_ignore = atoi(argv[3]);
    chains_to_use = num_chains - chains_ignore;

    Chain *chain_pdb;
    try {
        chain_pdb = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }


    for ( int i = 0; i < num_chains; i ++ ) {
        debug("building chain --> %d\n",i);
        strcpy(chain_pdb[i].filename,argv[1]);
        chain_pdb[i].chainid = i;
        if ( i == 0 ) {
            debug("building first chain (%d)\n", i);
            chain_pdb[i].index = 0;
            chain_pdb[i].file_line_begin = 0;
        } else {
            debug(" building later chain --> %d\n",i);
            chain_pdb[i].index = chain_pdb[i-1].findex + 1;
            chain_pdb[i].file_line_begin = chain_pdb[i-1].file_line_end; // or +1?
        }

        // Read PDB File.
        ReadEveryLine(&chain_pdb[i]);
        chain_pdb[i].findex = chain_pdb[i].index + chain_pdb[i].num_atoms - 1;
        chain_pdb[i].assign_indices();
        // chain_pdb[i].print_prop();

#ifdef INFO
        if (chain_ref[i].num_atoms_ca > 0) {chain_ref[i].print_prop();}
#endif // INFO
    }
    std::cout << "\n<<<<<<<< ALLATOM structure has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;
    std::cout << "Perform the memory reduction from PDB to AA (all atom).\n" << std::endl;

    // Coarse Grain. to CA
    Chain *chain_aa;
    try {
        chain_aa = new Chain [num_chains];
    } catch (std::bad_alloc xa) {
        std::cout << "Allocation Failure\n";
        return 1;
    }
    // exit(0);

    for ( int i = 0; i < num_chains; i ++ ) {
        // Undefine.
        debug("Redefining molecule size: (%d)\n",chain_pdb[i].num_atoms);

        // RIGHT SIZE !!!  <-- Read again. (know num_atoms).
#undef MOLECULE_INITIAL_SIZE
#define MOLECULE_INITIAL_SIZE chain_pdb[i].num_atoms

        debug("__building chain --> %d\n",i);
        strcpy(chain_aa[i].filename,argv[1]);
        chain_aa[i].chainid = i;
        if ( i == 0 ) {
            debug("__building first chain (%d)\n", i);
            chain_aa[i].index = 0;
            chain_aa[i].file_line_begin = 0;
        } else {
            debug("__building more chains --> %d\n",i);
            chain_aa[i].index = chain_aa[i-1].findex + 1;
            chain_aa[i].file_line_begin = chain_aa[i-1].file_line_end; // or +1?
        }
        // Read PDB File.
        ReadEveryLine(&chain_aa[i]);
        chain_aa[i].findex = chain_aa[i].index + chain_aa[i].num_atoms - 1;
        chain_aa[i].assign_indices();
#ifdef INFO
        if (chain_aa[i].num_atoms > 0) {chain_aa[i].print_prop();}
#endif // INFO
    }
    std::cout << " <<<<<<<< chain_aa has been acquired. # of chains: " << num_chains << " >>>>>>>> \n" << std::endl;
    // End of REFERENCE chain_aa


    // free
    // std::cout << " exiting safely .. sort of.\n" << std::endl;
    delete [] chain_pdb;

    // exit(0);
#endif // ALLATOM



    /* ---------------------------------------------------------
       dcd Preface complete.
       default args accepted.
       dcd for loop begun.
       --------------------------------------------------------- */

    /* ---------------------------------------------------------
       Processing.
       Post dcd loading.
       dcdloop
       All of the main functionality is described here.
       --------------------------------------------------------- */
    // * TODO, not complete
    // + attempted, not working.
    // _ done, working
    // D: define NAME
    // ---------------------------------------------------------
    // 1. Contact - generate, intra chain
    // 2. Contacts - by reference comparison
    // 3. Tension (in protein, 1-5 chains)
    // 4. Bond Vector, Angle with Tension   <-- same as 3.
    // 5.*chi:structural overlap
    // 6. Topology
    // 7. Angles (proto - only)
    // 8. Centroid of chain
    // 9. Centroid of centroids
    // 10. Microtubule contacts.
    // 11. Indices (get_or_print)


    // Check.
    // Vector m0;
    // Vector m1;
    // double dist;
    // m0.x = 94.574, m0.y = 6.982, m0.z = 126.656;
    // m1.x = 38.634, m1.y = 20.234,m1.z = 64.772;
    // dist = distance(m0,m1);
    // printf("the distance is: %f\n",dist);


#ifdef INDICES
    fprintf_all_indices(chain_ref,num_chains);
    // for ( int i=0; i<num_chains; i++ ) {
    //     // debug("chain:(%d) ",i);
    //     // build_contact_map(&chain_ref[i]);
    //     chain_ref[i].fprintf_indices();
    // }
#endif // INIDICES
#ifdef INDICES_AA
    fprintf_all_indices(chain_aa,num_chains);
#endif // INDICES_AA

#ifdef TOP_CHARMM

    // FILE
    FILE * fp_topology_charmm;
    fp_topology_charmm = fopen("topology_charmm.top", "w+");
    // fprintf(fp_topology_charmm,"\n");
    fprintf(fp_topology_charmm,"# new charmm topology:\n");

    top_charmm_write_general(chain_aa,num_chains,fp_topology_charmm);

    fclose(fp_topology_charmm);

#endif // TOP_CHARMM

#ifdef INERTIA
    // get Moment of Inertia Tensor
    // debug("")
    printf("entering inertia!\n");
    for ( int i=0; i<chains_to_use; i++) {
        debug("inertia: %d\n",i);
        chain_ref[i].GetMomentofInertiaTensor();
    }


#endif // INERTIA

#ifdef CURVATURE_M
    // std::cout << "begin CURVATURE evaluation." << std::endl;
    debug("begin CURVATURE evaluation.\n");

    // print Centroids
    for ( int i=0; i <chains_to_use; i++) {
        chain_later[i].ComputeCentroid();
    }

    // CRITERION: broken protofilament?
    // temporary > 50.0
    // permanent > 75.0

    double d1,d_max;
    d1 = d_max = 0.0;

    // start frame
    // fprintf(fp_proto_centroid_dist,"# start_frame: %d",frame_position);

    // only check the inter-dimer
    for(int i=1; i<chains_to_use-1; i+=2){
        d1 = distance(chain_later[i].centroid,chain_later[i+1].centroid);
        debug("distance_proto_centroids: (%d-%d) :>  %f\n",i,i+1,d1);
        fprintf(fp_proto_centroid_dist,"%7.2f  ",d1);
        if(d1 > d_max){
            d_max = d1;
        }
    }
    fprintf(fp_proto_centroid_dist,"\n");

    // printf("d_max: %f\n",d_max);
    // if ( d_max > 90.0){
    //     // FILE
    //     FILE * fp_framecount;
    //     // fp_framecount = fopen("framecount_%d.dat",(int)frame_position, "a+");
    //     fp_framecount = fopen("framecount_.dat","a+");
    //     // fprintf("%d\n",frame_position);
    //     fprintf(fp_framecount,"\n");
    //     fclose(fp_framecount);
    //     printf("leaving now: %d",frame_position);
    //     break;
    //     // sleep(1);
    //     // std::cin.ignore();
    // }

    // continue;
    // CRITERION: END.


    // chain system, total # chains,
    //   from axis, to axis,begin point <-> end point (vector)
    // void rotate_system_around_axis(Chain *chain,int max_num_chains,
    //                                        int axisfrom,int axisto,Vector O,Vector E);


    /* ---------------------------------------------------------
       Traditional Rotation of entire system.
       --------------------------------------------------------- */
    Vector t1;
    Vector newpos;

    // translate to origin
    t1 = translate_to_origin(chain_later,0,num_chains);
    /* ---------------------------------------------------------
       Rotation of the entire system completed.
       --------------------------------------------------------- */

    // Begin global rotating._____________________________________________
    // z_rotate
    for(int i=0; i<num_chains; i++){
        for(int j=0; j<chain_later[i].num_atoms_ca; j++){
            // printf("%d",j);
            newpos = prod_matrix_vector(Rz,chain_later[i].pos[j]);
            chain_later[i].pos[j].x = newpos.x;
            chain_later[i].pos[j].y = newpos.y;
            chain_later[i].pos[j].z = newpos.z;
        }
    }
    // x_rotate
    for(int i=0; i<num_chains; i++){
        for(int j=0; j<chain_later[i].num_atoms_ca; j++){
            // printf("%d",j);
            newpos = prod_matrix_vector(Rx,chain_later[i].pos[j]);
            chain_later[i].pos[j].x = newpos.x;
            chain_later[i].pos[j].y = newpos.y;
            chain_later[i].pos[j].z = newpos.z;
        }
    }
    // y_rotate
    for(int i=0; i<num_chains; i++){
        for(int j=0; j<chain_later[i].num_atoms_ca; j++){
            // printf("%d",j);
            newpos = prod_matrix_vector(Ry,chain_later[i].pos[j]);
            chain_later[i].pos[j].x = newpos.x;
            chain_later[i].pos[j].y = newpos.y;
            chain_later[i].pos[j].z = newpos.z;
        }
    }
    // rotated about Y, from yz quad 1, to xy, quad 1.
    Matrix Rot;
    Rot.build_rotation(M_PI * 0.5 * -1.0,1); // theta,Y:1 ... -90 d about y.
    // Rot.print_Matrix();

    // rotate
    // for(int i=0; i<chains_to_use; i++){
    for(int i=0; i<num_chains; i++){
        for(int j=0; j<chain_later[i].num_atoms_ca; j++){
            // printf("%d",j);
            newpos = prod_matrix_vector(Rot,chain_later[i].pos[j]);
            chain_later[i].pos[j].x = newpos.x;
            chain_later[i].pos[j].y = newpos.y;
            chain_later[i].pos[j].z = newpos.z;
        }
    }
    // End of global rotating.____________________________________________



    // get global curvature
    Vector topcurve_midpoint;
    int mid_chain1,mid_chain2;
    mid_chain2 = chains_to_use * 0.5;
    mid_chain1 = mid_chain2 - 1;
    // printf("mid_chain1: %d %d\n",mid_chain1,mid_chain2);
    // get_vector(chain_later[mid_chain1].centroid,chain_later[mid_chain2].centroid,&topcurve_midpoint);
    topcurve_midpoint = midpoint2(chain_later[mid_chain1].centroid,chain_later[mid_chain2].centroid);


    double dt1,dt2;
    dt1 = distance(chain_later[mid_chain1].centroid,topcurve_midpoint);
    dt2 = distance(chain_later[mid_chain2].centroid,topcurve_midpoint);
    debug("distances(to p2): either  %f  OR   %f\n",dt1,dt2);

    // degress, curvature variable;
    Curvature global_curv;
    Curvature local_curv;

    // GLOBAL: get global curvature
    debug("\nglobal_curvature:\n");
    global_curv.p1 = chain_later[0].centroid;
    global_curv.p2 = topcurve_midpoint;
    global_curv.p3 = chain_later[chains_to_use-1].centroid;


#ifndef NDEBUG
    //
    if ((global_curv.radius != 0.0) || (global_curv.curvature != 0.0)){
        printf("Warning! curvature & radius aren't zero.\n");
        global_curv.print_points();
        global_curv.print_rc();
        exit(1);
    }
#endif

#ifdef CURV4
    global_curv = compute_curvature4(chain_later,0,mid_chain1,mid_chain2,chains_to_use-1,fp_curvature_global,fp_radius_of_curvature_global);

    // NEW
    fprintf(fp_curvature_global,"\n");
    fprintf(fp_radius_of_curvature_global,"\n");

    for(int i=0; i<=chains_to_use-4; i+=2) {
        printf("\nlocal_curvature: %d\n",i);
        // local_curv = compute_curvature(chain_later,i,fp_curvature_local,fp_radius_of_curvature_local);
        local_curv = compute_curvature4(chain_later,i,i+1,i+2,i+3,fp_curvature_local,fp_radius_of_curvature_local);
        // local_curv.print_points(); // avgcurv returned, no points.
        // local_curv.print_rc();
    }
    // FPRINTF - local.
    fprintf(fp_curvature_local,"\n");
    fprintf(fp_radius_of_curvature_local,"\n");
#endif // CURV4



#ifdef CURV3
    global_curv = compute_curvature3(chain_later,0,mid_chain2,chains_to_use-1,fp_curvature_global,fp_radius_of_curvature_global);

    // NEW
    fprintf(fp_curvature_global,"\n");
    fprintf(fp_radius_of_curvature_global,"\n");

    for(int i=0; i<=chains_to_use-2; i+=2) {
        printf("\nlocal_curvature: %d\n",i);
        // local_curv = compute_curvature(chain_later,i,fp_curvature_local,fp_radius_of_curvature_local);
        local_curv = compute_curvature3(chain_later,i,i+1,i+2,fp_curvature_local,fp_radius_of_curvature_local);
        // local_curv.print_points(); // avgcurv returned, no points.
        // local_curv.print_rc();
    }
    // FPRINTF - local.
    fprintf(fp_curvature_local,"\n");
    fprintf(fp_radius_of_curvature_local,"\n");
#endif // CURV3


    // debug("<-- quality_check -->\n");
    // Vector p1,p2,p3;
    // p1.x = 1.0;
    // p1.y = 1.0;
    // p2.x = 2.0;
    // p2.y = 3.0;
    // p3.x = 3.0;
    // p3.y = 8.0;
    // curvature(p1,p2,p3);
    // exit(0);
#endif // CURVATURE_M


#ifdef CONTACTMAP_SM2
    fprintf(fp_contacts,"# frame_position: %d\n",frame_position);


#ifdef CONTACTSBYRESFILE1
    // fp_contactsresfile1 = fopen("contacts_byresidue.dat");
    fprintf(fp_contactsresfile1,"# frame_position: %d\n",frame_position);
    topo_contacts_persisting_by_residue(chain_later,0,chain_ref[0].contacts,fp_contactsresfile1);
#endif // CONTACTSBYRESFILE1


    int contacts = 0;
    // for (int i=0; i<chains_to_use; i++) {
    topo_count_map(chain_ref[0].contacts);
    // contacts = topo_contacts_persisting(chain_later,0,chain_ref[0].contacts,20.0);
    contacts = topo_contacts_persisting(chain_later,0,chain_ref[0].contacts,fp_contacts);
    // fprintf(fp_contacts," %d",contacts);

    // }
    // FILE * fp_contacts;
    // fp_contacts_n = fopen("contacts_intra.dat", "w+");




#endif // CONTACTMAP_SM2

#ifdef CONTACTMAP_MT2

    // N
    // Contact map_neighbor[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
    // for(int i=0;i<chains_to_use; i++){
    //     // topo_build_inter(chain_later,i,chain_ref[i].LonN,map_neighbor,20.0);
    //     memcpy(chain_later[i].contactsLonN,map_neighbor, sizeof(chain_later[i].contactsLonN));
    //     for(int j=0; j<MAX_CONTACTS; j++){
    //         map_neighbor[j]->~Contact();
    //     }
    // }

    int contacts = 0;

#ifdef CONFILE
    fprintf(fp_contacts_n,"# frame: %d\n",frame_position);
    fprintf(fp_contacts_e,"# frame: %d\n",frame_position);
    fprintf(fp_contacts_s,"# frame: %d\n",frame_position);
    fprintf(fp_contacts_w,"# frame: %d\n",frame_position);
#endif // CONMAP


    std::vector< std::vector <int> > persisting_contacts;
    std::vector<int> persist_nesw;

    // for (it=pf_south.begin(); it<pf_south.end(); it++) {


    // new lines.
    // if(chains_to_use > 1){
    for (int i=0; i<chains_to_use; i++) {

        // N (the 20.0 was when it was a soft criterion)
        // topo_count_map(chain_ref[i].contactsLonN);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LonN,chain_ref[i].contactsLonN,20.0);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LonN,chain_ref[i].contactsLonN);
        contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LonN,chain_ref[i].contactsLonN,13.0); // cutoff
        persist_nesw.insert(persist_nesw.end(),contacts);
#ifdef DIMERMAP
        maprow_contacts.push_back(contacts);
#endif
#ifdef CONFILE
        fprintf(fp_contacts_n," %d",contacts);
#endif // CONFILE

        // E
        // topo_count_map(chain_ref[i].contactsLatE);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LatE,chain_ref[i].contactsLatE,20.0);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LatE,chain_ref[i].contactsLatE);
        contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LatE,chain_ref[i].contactsLatE,13.0); // cutoff
        persist_nesw.insert(persist_nesw.end(),contacts);
#ifdef DIMERMAP
        maprow_contacts.push_back(contacts);
#endif
#ifdef CONFILE
        fprintf(fp_contacts_e," %d",contacts);
#endif // CONFILE


        // S
        // topo_count_map(chain_ref[i].contactsLonS);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LonS,chain_ref[i].contactsLonS,20.0);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LonS,chain_ref[i].contactsLonS);
        contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LonS,chain_ref[i].contactsLonS,13.0); // cutoff
        persist_nesw.insert(persist_nesw.end(),contacts);
#ifdef DIMERMAP
        maprow_contacts.push_back(contacts);
#endif
#ifdef CONFILE
        fprintf(fp_contacts_s," %d",contacts);
#endif // CONFILE


        // W
        // topo_count_map(chain_ref[i].contactsLatW);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LatW,chain_ref[i].contactsLatW,20.0);
        // contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LatW,chain_ref[i].contactsLatW);
        contacts = topo_contacts_persisting(chain_later,i,chain_ref[i].LatW,chain_ref[i].contactsLatW,13.0);
        persist_nesw.insert(persist_nesw.end(),contacts);
#ifdef DIMERMAP
        maprow_contacts.push_back(contacts);
#endif
#ifdef CONFILE
        fprintf(fp_contacts_w," %d",contacts);
#endif /// CONFILE
        // count_map(map_ref);
        // clear_map(map_ref);


        persisting_contacts.push_back(persist_nesw);
        persist_nesw.clear();


#ifdef DIMERMAP
        map_contacts.insert( std::pair<int, std::vector<int> >(i,maprow_contacts));
        maprow_contacts.clear();
#endif
    }

    // fprintf(fp_contacts_n,"\n# frame: %d\n",frame_position);
    // fprintf(fp_contacts_e,"\n# frame: %d\n",frame_position);
    // fprintf(fp_contacts_s,"\n# frame: %d\n",frame_position);
    // fprintf(fp_contacts_w,"\n# frame: %d\n",frame_position);

#ifdef CONFILE
    fprintf(fp_contacts_n,"\n");
    fprintf(fp_contacts_e,"\n");
    fprintf(fp_contacts_s,"\n");
    fprintf(fp_contacts_w,"\n");
#endif // CONFILE
    // }


    // std::vector<int>::const_iterator it3;
    // std::vector<int>::iterator it2;
    // for(int it3=0; it3<persisting_contacts.size(); it3++){
    //     // printf("%d ",it3);
    //     // for(int it4=0; it4<persisting_contacts
    //     printf("%d %d %d %d\n",persisting_contacts[0],\
    //            persisting_contacts[1], \
    //            persisting_contacts[2], \
    //            persisting_contacts[3]);
    // }


// #ifndef NDEBUG
//     // check the 2D array: pf_array.
//     printf("2d_array:\n");
//     for (int i=0; i<persisting_contacts.size(); i++) {
//         printf("%d:\t",i);
//         for (int j=0; j<persisting_contacts[i].size(); j++) {
//             printf("%3d ",persisting_contacts[i][j]);
//         }
//         printf("\n");
//     }
// #endif


    fprintf(fp_dimercontacts,"# frame: %d\n",frame_position);

    // DCD: assemble dimer contact counts.
    for (int numdimers=0; numdimers<dimers_alpha.size(); numdimers++) {
        debug("%d ---> %d\n",numdimers,dimers_alpha[numdimers]);


        topo_sort_dimer_contacts_initially(chain_later,chain_ref,   \
                                           dimers_alpha[numdimers], \
                                           fp_dimercontacts,        \
                                           persisting_contacts);

    }
#endif // CONTACTMAP_MT2

#ifdef MTPF_M
    printf("MTPF_M: processing..\n");

    int ch;
    int contacts = 0;

    for ( int k = 0; k < pf_array.size()*4; k++ ) {
        outcontactfile[k] << "# frame: " << frame_position << "\n";
    }

    // evaluate contacts.
    for (int i = 0; i<pf_array.size(); i++) {
        debug("getting contacts:\t%d\n",i);

        for (int d = 0; d<4; d++ ) {

            for (int j = 0; j<pf_array[i].size(); j++) {

                ch = pf_array[i][j]; //
                printf("%3d ",ch);

                // N,E,S,W
                if ( d == 0 ) {
                    contacts = topo_contacts_persisting(chain_later,ch,chain_ref[ch].LonN,chain_ref[ch].contactsLonN);
                } else if ( d == 1 ) {
                    contacts = topo_contacts_persisting(chain_later,ch,chain_ref[ch].LatE,chain_ref[ch].contactsLatE);
                } else if ( d == 2 ) {
                    contacts = topo_contacts_persisting(chain_later,ch,chain_ref[ch].LonS,chain_ref[ch].contactsLonS);
                } else {
                    contacts = topo_contacts_persisting(chain_later,ch,chain_ref[ch].LatW,chain_ref[ch].contactsLatW);
                }

                outcontactfile[i*4+d] << std::setw(3) << contacts << ' ';

            } // j monomer members

            outcontactfile[i*4+d] << '\n';

        } // d: 0-3, for the N E S W

        printf("\n");
    } // i: 0-12 if 13 protofilaments.

    printf("MTPF_M: complete.\n");
#endif // MTPF_M

#ifdef MT_DIMER_CONTACT_M


#endif // MT_DIMER_CONTACT_M

#ifdef MTTOP
    // std::cout << "MT! using chains: " << chains_to_use << std::endl;
    debug("MTTOP: chains: %d\n",chains_to_use);

    Vector ref_centroid_pos[chains_to_use-1];
    // std::vector<Vector> ref_centroid_pos[chains_to_use-1];

    for ( int i=0; i<chains_to_use; i ++ ) {
        get_tubulin_centroid(&chain_ref[i],&ref_centroid_pos[i]);
        // get_tubulin_centroid(&chain_ref[i],&ref_centroid_pos[i]);
    }
    Vector centroid_front;
    Vector centroid_back;
    Vector mt_axis;
    Vector mt_axis_norm;
    // centroid_front.x = 0.0, centroid_front.y = 0.0, centroid_front.z = 0.0;
    // std::cout << centroid_front.x << centroid_front.y << centroid_front.z << std::endl;
    // std::cout << centroid_back.x << centroid_back.y << centroid_back.z << std::endl;
    // std::cout << mt_axis.x << mt_axis.y << mt_axis.z << std::endl;
    // std::cout << mt_axis_norm.x << mt_axis_norm.y << mt_axis_norm.z << std::endl;

    get_centroid_from_array(ref_centroid_pos, 0, chains_to_use/2, &centroid_front);
    // get_centroid_from_array(ref_centroid_pos, 0, 30, &centroid_front);
    get_centroid_from_array(ref_centroid_pos, chains_to_use/2 + 1, chains_to_use, &centroid_back);
    // get_centroid_from_array(ref_centroid_pos, chains_to_use - 30, chains_to_use, &centroid_back);
    get_vector(centroid_front, centroid_back, &mt_axis);
    normalize(mt_axis, &mt_axis_norm);
    debug("Centroid_front: %f %f %f\n",centroid_front.x,centroid_front.y,centroid_front.z);
    debug("Centroid_back: %f %f %f\n",centroid_back.x,centroid_back.y,centroid_back.z);
    debug("mt_axis: %f %f %f\n",mt_axis.x,mt_axis.y,mt_axis.z);
    debug("mt_axis_norm: %f %f %f\n",mt_axis_norm.x,mt_axis_norm.y,mt_axis_norm.z);


    // rewrite the centroid, to use chain_ref[i].centroid.x.y.z;
    printf("new way to get centroids ..\n");
    for ( int i=0; i < num_chains; i++) {
        chain_ref[i].ComputeCentroid();
    }

    printf(".. extending analysis to contact interfaces .......\n");
    compute_nearest_neighbors(chain_ref,chains_to_use);

    printf(".. determining latitudinal and longitudinal neighbors .......\n");
    determine_latlon_neighbors(chain_ref,chains_to_use,mt_axis_norm);

    for ( int i=0; i < num_chains; i++) {
    // for ( int i=5; i < 15; i++) {
        debug("chain:(%d) \n",i);
        debug("\n%3d >> N:%3d S:%3d W:%3d E:%3d\n",i,chain_ref[i].LonN,chain_ref[i].LonS,chain_ref[i].LatW,chain_ref[i].LatE);
        printf("chain:(%d)\n",i);
        printf("%3d >>N:%3d E:%3d S:%3d W:%3d\n",i,chain_ref[i].LonN,chain_ref[i].LatE,chain_ref[i].LonS,chain_ref[i].LatW);
    }

    printf("building intra __contacts__ map..\n");
    for ( int i=0; i<chains_to_use; i++ ) {
        debug("chain:(%d) \n",i);
        build_contact_map(&chain_ref[i]);
    }

    printf(".. getting latitudinal and longitudinal contacts .......\n");
    // N E S W | 1 2 3 4
    for ( int i=0; i<chains_to_use; i++ ) {

        Contact map1[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter(chain_ref,i,chain_ref[i].LonN,map1); // important: 1,0 small/big
        memcpy(chain_ref[i].contactsLonN, map1, sizeof(chain_ref[i].contactsLonN));
        Contact map2[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter(chain_ref,i,chain_ref[i].LatE,map2); // important: 1,0 small/big
        memcpy(chain_ref[i].contactsLatE, map2, sizeof(chain_ref[i].contactsLatE));
        Contact map3[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter(chain_ref,i,chain_ref[i].LonS,map3); // important: 1,0 small/big
        memcpy(chain_ref[i].contactsLonS, map3, sizeof(chain_ref[i].contactsLonS));
        Contact map4[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter(chain_ref,i,chain_ref[i].LatW,map4); // important: 1,0 small/big
        memcpy(chain_ref[i].contactsLatW, map4, sizeof(chain_ref[i].contactsLatW));

        // memcpy(chain_ref[i].contactsLonN, map1, sizeof(chain_ref[i].contactsLonN));
        // Contact map2[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // build_contact_map_inter(chain_ref,i,map2,2);
        // memcpy(chain_ref[i].contactsLatE, map2, sizeof(chain_ref[i].contactsLatE));
        // Contact map3[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // build_contact_map_inter(chain_ref,i,map3,3);
        // memcpy(chain_ref[i].contactsLonS, map3, sizeof(chain_ref[i].contactsLonS));
        // Contact map4[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // build_contact_map_inter(chain_ref,i,map4,4);
        // memcpy(chain_ref[i].contactsLatW, map4, sizeof(chain_ref[i].contactsLatW));
    }


    printf("assigning eh..\n");
    for ( int i=0; i<chains_to_use; i++ ) {
        debug("chain:(%d) \n",i);
        assign_contact_eh(chain_ref[i].contacts,1.33333); // gets the rest of the 0.0;
        assign_contact_eh(chain_ref[i].contactsLonN,1.99999);
        assign_contact_eh(chain_ref[i].contactsLatE,1.88888);
        assign_contact_eh(chain_ref[i].contactsLonS,1.99999);
        assign_contact_eh(chain_ref[i].contactsLatW,1.88888);
    }

    printf("printing topology..\n");
    for ( int i=0; i<chains_to_use; i++ ) {
    // for ( int i=0; i<5; i++ ) {

        printf("(%d) I:      intra ",i);
        fprintf_contact_map_gsop(chain_ref[i].contacts);

        if ( i < chain_ref[i].LonN) {
            printf("N: ");
            fprintf_contact_map_gsop(chain_ref[i].contactsLonN);
        }
        if ( i < chain_ref[i].LatE) {
            printf("E: ");
            fprintf_contact_map_gsop(chain_ref[i].contactsLatE);
        }
        if ( i < chain_ref[i].LonS) {
            printf("S: ");
            fprintf_contact_map_gsop(chain_ref[i].contactsLonS);
        }
        if ( i < chain_ref[i].LatW) {
            printf("W: ");
            fprintf_contact_map_gsop(chain_ref[i].contactsLatW);
        }
    }

    for (int i=0; i<chains_to_use; i++)
    {
        printf("N:%3d  E:%3d  S:%3d  W:%3d\n",chain_ref[i].LonN_total2x * 0.5,\
               chain_ref[i].LatE_total2x * 0.5,                         \
               chain_ref[i].LonS_total2x * 0.5,                         \
               chain_ref[i].LatW_total2x * 0.5);
    }


#endif // MTTOP

#ifdef TOPOPF
    printf("\nhello from protofilament topology writer!\n\n");

    Contact map[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
    // Contact map_interN[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
    // Contact map_interE[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
    // Contact map_interS[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
    // Contact map_interW[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];


    // intra --> contacts
    for(int i=0; i<num_chains; i++)
    {
        chain_ref[i].ComputeCentroid();
        topo_count_map(map);
        topo_build_intra(chain_ref,i,map,8.0); // *chain,int chain_num,Contact (*map)[MAX_CONTACTS],cutoff
        memcpy(chain_ref[i].contacts,map,sizeof(chain_ref[i].contacts));
        topo_count_map(map);
        topo_clear_map(map);
        assign_contact_eh(chain_ref[i].contacts,2.00);
    }

    // Southern --> map --> contactsLonS
    for(int i=1; i<chains_to_use; i++)
    {
        topo_count_map(map);
        topo_build_inter(chain_ref,i,i-1,map,8.0);
        memcpy(chain_ref[i].contactsLonS,map,sizeof(chain_ref[i].contactsLonS));
        topo_count_map(map);
        topo_clear_map(map);
        assign_contact_eh(chain_ref[i].contactsLonS,2.00);
    }

    // Southern --> map --> contactsLonS
    for(int i=0; i<chains_to_use-1; i++)
    {
        topo_count_map(map);
        topo_build_inter(chain_ref,i,i+1,map,8.0);
        memcpy(chain_ref[i].contactsLonN,map,sizeof(chain_ref[i].contactsLonN));
        topo_count_map(map);
        topo_clear_map(map);
        assign_contact_eh(chain_ref[i].contactsLonN,2.00);
    }

    // i kinesins.
    float kd,kdmin;
    kd = kdmin = 100.0;
    int neighbor1,neighbor2;
    neighbor1 = neighbor2 = -1;

    for(int i=chains_to_use; i<num_chains; i++)
    {
        printf("kinesin: %d\n",i);
        kd = kdmin = 100.0;
        neighbor1 = neighbor2 = -1;

        // j alpha/beta tubulin.
        for(int j=0; j<chains_to_use; j++)
        {
            kd = distance (chain_ref[i].centroid,chain_ref[j].centroid);
            // printf("the distance between kinesin and the tubulin was: %f (%d)\n",kd,j);
            if(kd < kdmin)
            {
                kdmin = kd;
                neighbor1 = j;
            }
        }


        kd = kdmin = 100.0;

        for(int j=0; j<chains_to_use; j++)
        {
            kd = distance (chain_ref[i].centroid,chain_ref[j].centroid);
            // printf("the distance between kinesin and the tubulin was: %f (%d)\n",kd,j);
            if((kd < kdmin) and (j != neighbor1))
            {
                kdmin = kd;
                neighbor2 = j;
            }
        }

        printf("the closest neighbor was: %d\n",neighbor1);
        printf("the second closest neighbor was: %d\n",neighbor2);

        // chain_ref[i].LatW = neighbor1;
        // chain_ref[i].LatE = neighbor2;

        topo_count_map(map);
        topo_build_inter(chain_ref,i,neighbor1,map,8.0);
        memcpy(chain_ref[i].contactsLatW,map,sizeof(chain_ref[i].contactsLatW));
        topo_count_map(map);
        topo_clear_map(map);
        assign_contact_eh(chain_ref[i].contactsLatW,2.00);

        topo_count_map(map);
        topo_build_inter(chain_ref,i,neighbor2,map,8.0);
        memcpy(chain_ref[i].contactsLatE,map,sizeof(chain_ref[i].contactsLatE));
        topo_count_map(map);
        topo_clear_map(map);
        assign_contact_eh(chain_ref[i].contactsLatE,2.00);
    }

    // exit(0);


    // Value assignment:
    for(int i=0; i<chains_to_use; i++)
    {
        fprintf_contact_map_gsop(chain_ref[i].contacts);
        fprintf_contact_map_gsop(chain_ref[i].contactsLonS);
        fprintf_contact_map_gsop(chain_ref[i].contactsLonN);
    }

    for(int i=chains_to_use; i<num_chains; i++)
    {
        fprintf_contact_map_gsop(chain_ref[i].contacts);
        fprintf_contact_map_gsop(chain_ref[i].contactsLatW); // for near neighbor 1.
        fprintf_contact_map_gsop(chain_ref[i].contactsLatE); // for near neighbor 2.
    }


    // print_contact_map()
    // print_contact_map(chain_ref[0].contacts); // 1176
    // print_contact_map(chain_ref[1].contacts); // happens to be 0.
    // print_contact_map(chain_ref[0].intercontacts); // 80
    // print_contact_map(chain_ref[1].intercontacts); // same.

    // print gsop topology, native section.
    // ---------------------------------------------------------------------


    // FILE * fp_contacts;
    // fp_contacts = fopen("contacts_intra.dat", "w+");

#endif // TOPOPF


#ifdef TOPOLOGY
    for ( int i=0; i < num_chains; i++ ) {
        debug("chain:(%d) \n",i);
        build_contact_map(&chain_ref[i]);
    }
    // BIG then SMALL failed. (ca:216 with ca:7)
    // Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    // // void build_contact_map_inter(Chain *chain_later,int cid1,int cid2,Contact (*map)[MAX_CONTACTS]);
    // build_contact_map_inter(chain_ref,0,1,map);
    // memcpy(chain_ref[0].intercontacts, map, sizeof(chain_ref[0].intercontacts));
    // ---------------------------------------------------------------------

#ifdef SBDPEP
    // SMALL then BIG worked.  (ca:7 with ca:216)
    Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    build_contact_map_inter(chain_ref,1,0,map); // important: 1,0 small/big
    memcpy(chain_ref[0].intercontacts, map, sizeof(chain_ref[0].intercontacts));
    // memcpy(chain_ref[1].intercontacts, map, sizeof(chain_ref[1].intercontacts));
    // ---------------------------------------------------------------------
    // print_contact_map()
    // print_contact_map(chain_ref[0].contacts); // 1176
    // print_contact_map(chain_ref[1].contacts); // happens to be 0.
    // print_contact_map(chain_ref[0].intercontacts); // 80
    // print_contact_map(chain_ref[1].intercontacts); // same.

    // assign eh.
    // assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh,int low_resid,\
    //                        int high_resid,int low_cresid,int high_cresid);
    // ---------------------------------------------------------------------
    assign_contact_eh(chain_ref[0].contacts,1.35,0,113,0,113);
    assign_contact_eh(chain_ref[0].contacts,1.80,114,215,114,215);
    // print_contact_map(chain_ref[0].contacts);
    // debug("contact mapping -------------");
    assign_contact_eh(chain_ref[0].contacts,1.55885); // gets the rest of the 0.0;
    // print_contact_map(chain_ref[0].contacts);
    // debug("contact mapping -------------");
    assign_contact_eh(chain_ref[0].intercontacts,1.55885);


    // print gsop topology, native section.
    // ---------------------------------------------------------------------
    fprintf_contact_map_gsop(chain_ref[0].contacts);
    fprintf_contact_map_gsop(chain_ref[1].contacts);
    fprintf_contact_map_gsop(chain_ref[0].intercontacts);
#endif // SBDPEP
#ifdef SBDNOPEP
    // SMALL then BIG worked.  (ca:7 with ca:216)
    // Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    // build_contact_map_inter(chain_ref,1,0,map); // important: 1,0 small/big
    // memcpy(chain_ref[0].intercontacts, map, sizeof(chain_ref[0].intercontacts));
    // memcpy(chain_ref[1].intercontacts, map, sizeof(chain_ref[1].intercontacts));
    // ---------------------------------------------------------------------
    // print_contact_map()
    // print_contact_map(chain_ref[0].contacts); // 1176
    // print_contact_map(chain_ref[1].contacts); // happens to be 0.
    // print_contact_map(chain_ref[0].intercontacts); // 80
    // print_contact_map(chain_ref[1].intercontacts); // same.

    // assign eh.
    // assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh,int low_resid,\
    //                        int high_resid,int low_cresid,int high_cresid);
    // ---------------------------------------------------------------------
    assign_contact_eh(chain_ref[0].contacts,1.35,0,119,0,119);
    assign_contact_eh(chain_ref[0].contacts,1.80,120,221,120,221);
    assign_contact_eh(chain_ref[0].contacts,1.55885); // gets the rest of the 0.0;
    // print_contact_map(chain_ref[0].contacts);
    // debug("contact mapping -------------");

    // print_contact_map(chain_ref[0].contacts);
    // debug("contact mapping -------------");
    // assign_contact_eh(chain_ref[0].intercontacts,1.55885);


    // print gsop topology, native section.
    // ---------------------------------------------------------------------
    fprintf_contact_map_gsop(chain_ref[0].contacts);
    // fprintf_contact_map_gsop(chain_ref[1].contacts);
    // fprintf_contact_map_gsop(chain_ref[0].intercontacts);
#endif // SBDNOPEP


// #ifdef MTTOP
//     // SMALL then BIG worked.  (ca:7 with ca:216)

//     for ( int i=0; i < num_chains; i++ ) {
//         debug("chain:(%d) ",i);
//         build_contact_map(&chain_ref[i]);
//     }

//     Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
//     build_contact_map_inter(chain_ref,1,0,map); // important: 1,0 small/big
//     memcpy(chain_ref[0].intercontacts, map, sizeof(chain_ref[0].intercontacts));
//     // memcpy(chain_ref[1].intercontacts, map, sizeof(chain_ref[1].intercontacts));
//     // ---------------------------------------------------------------------
//     // print_contact_map()
//     // print_contact_map(chain_ref[0].contacts); // 1176
//     // print_contact_map(chain_ref[1].contacts); // happens to be 0.
//     // print_contact_map(chain_ref[0].intercontacts); // 80
//     // print_contact_map(chain_ref[1].intercontacts); // same.

//     // assign eh.
//     // assign_contact_eh(Contact (*map)[MAX_CONTACTS],float eh,int low_resid,\
//     //                        int high_resid,int low_cresid,int high_cresid);
//     // ---------------------------------------------------------------------
//     assign_contact_eh(chain_ref[0].contacts,1.35,0,113,0,113);
//     assign_contact_eh(chain_ref[0].contacts,1.80,114,215,114,215);
//     // print_contact_map(chain_ref[0].contacts);
//     // debug("contact mapping -------------");
//     assign_contact_eh(chain_ref[0].contacts,1.55885); // gets the rest of the 0.0;
//     // print_contact_map(chain_ref[0].contacts);
//     // debug("contact mapping -------------");
//     assign_contact_eh(chain_ref[0].intercontacts,1.55885);


//     // print gsop topology, native section.
//     // ---------------------------------------------------------------------
//     fprintf_contact_map_gsop(chain_ref[0].contacts);
//     fprintf_contact_map_gsop(chain_ref[1].contacts);
//     fprintf_contact_map_gsop(chain_ref[0].intercontacts);
// #endif // MTTOP

#endif // TOPOLOGY



// #ifdef CONTACT
//     // Section 1. Contacts (intra only)
//     for ( int i=0; i < num_chains; i++) {
//         printf("chain:(%d) ",i);
//         build_contact_map(&chain_ref[i]);
//     }
//     for ( int i=0; i < num_chains; i++ ) {
//         printf("chain:(%d) ",i);
//         build_contact_map(&chain_later[i]);
//     }

// #ifdef CONTACTPERSIST
//     fprintf_count_of_contacts_2x_persist_from_reference(&chain_ref[0],&chain_later[0]);
// #endif // CONTACTPERSIST
// #endif // CONTACT



    /* ---------------------------------------------------------
       Contact map comparison preparation
       bond_vector_angle_with_tension_vector (costheta, tension)
       --------------------------------------------------------- */
#ifdef TENSION_COSTHETA
    fprintf_bond_vector_angle_with_tension_vector(&chain_ref[0],&chain_later[0]);
#endif // TENSION_COSTHETA

    /* ---------------------------------------------------------
       Chi Value.
       --------------------------------------------------------- */
#ifdef CHI_MID
    // get_chi(chain_ref,chain_later,181,540);
    double chi;
    chi = 0.0;

    // chi = get_chi_global(chain_ref,chain_later,limit_begin,limit_end);
    // fprintf(fp_chi,"%4.3f\n",chi);

    // fprint(fp_chi,"# frame:"
    fprintf(fp_chi,"# frame: %d\n",frame_position);
    get_chi_by_residue(chain_ref,chain_later,limit_begin,limit_end,fp_chi);


#endif // Chi Value.



    /* ---------------------------------------------------------
       Protofilament Angle Analysis
       --------------------------------------------------------- */
#ifdef PFANGLE
    // usage
    if (argc != 8) {
        std::cout << "Usage: " << argv[0] \
                  << " <Filename-PDB>" \
                  << " <Filename-DCD>" \
                  << " <step (for DCD)>" \
                  << " <num_chains>" \
                  << " <ignore_chains>" \
                  << " <atom1>" \
                  << " <atom2>" \
                  << "\ncreates: angle_tubulins_byresids.dat (a+)" \
                  << std::endl;
        exit(1);
    }
    int resid1 = atoi(argv[6]);
    int resid2 = atoi(argv[7]);



    /* evaluate angles */
    FILE * fp_tubulin_monomers_angles;
    fp_tubulin_monomers_angles = fopen("angle_tubulins_byresids.dat", "a+");
    fprintf(fp_tubulin_monomers_angles,"# frame: %d\n",frame_position);

    for ( int i=0; i < num_chains-1-chains_ignore; i++) {
        /* -3 because it's minus 1, and not 2 kinesins */
        debug("\n%d %d\n",chain_later[i].num_atoms_ca,chain_later[i+1].num_atoms_ca);
        tubulin_monomers_angles(&chain_later[i],&chain_later[i+1],resid1,resid2,fp_tubulin_monomers_angles);
    }

    fprintf(fp_tubulin_monomers_angles,"\n");
    fclose(fp_tubulin_monomers_angles);
#endif // PFANGLE


#ifdef CENTROIDMOVEMENT_M
    printf("chain_ref map of networked/connected monomers: set map for chain_later, get contacts......\n");
    for ( int i=0; i<chains_to_use; i++ ) {
        // N E S W | 1 2 3 4
        // printf(">> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);
        chain_later[i].LonN = chain_ref[i].LonN;
        chain_later[i].LatE = chain_ref[i].LatE;
        chain_later[i].LonS = chain_ref[i].LonS;
        chain_later[i].LatW = chain_ref[i].LatW;
        // printf("R> %d %d %d %d\n",chain_ref[i].LonN,chain_ref[i].LatE,chain_ref[i].LonS,chain_ref[i].LatW);
        // printf("L> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);
    }

    // ComputeCentroid()
    for ( int i=0; i < num_chains; i++) {
        chain_later[i].ComputeCentroid();
    }

    // compare centroid to centroid
    get_centroid_movement(chain_ref,chain_later,num_chains);
    get_centroid_xyzposition_movement(chain_later,chains_to_use,fp_xyzpos);
#endif // CENTROIDMOVEMENT_M



#ifdef ANGLE3CENTROID_M
    debug("entering ANGLE3CENTROID!!!!\n");

    // failed!! use chain_ref!!
    // for( int i=0; i<chains_to_use; i++) {
    //     chain_later[i].LonN = chain_ref[i].LonN;
    //     chain_later[i].LatE = chain_ref[i].LatE;
    //     chain_later[i].LonS = chain_ref[i].LonS;
    //     chain_later[i].LatW = chain_ref[i].LatW;
    // }

    // debug("\n\ngetting latangle from chain_later!\n\n");

    // FILE * fp_angles;
    // fp_angles_ns = fopen ("mt_angles_ns.dat", "a+");
    // debug("opening for # time later, mt_angles_ns.dat\n");
    // fprintf(fp_angles_ns, "# time later Center|West|East or C|South|North (frame:%d)\n",frame_position);
    // fclose(fp_angles_ns);

    // fp_angles_ew = fopen ("mt_angles_ew.dat", "a+");
    // debug("opening for # time later, mt_angles_ew.dat\n");
    // fprintf(fp_angles_ew, "# time later Center|West|East or C|South|North (frame:%d)\n",frame_position);
    // fclose(fp_angles_ew);



    for( int i=0; i<chains_to_use; i++) {
        chain_later[i].ComputeCentroid();
    }

    for( int i=0; i<chains_to_use; i++) {
        // debug("entering .. %d\n",i);
        get_latangle_from35centroids(chain_ref,chain_later,i,fp_angles_ns,0); // S,N 1-5
        get_latangle_from35centroids(chain_ref,chain_later,i,fp_angles_ew,1); // W,E 1-5
        // debug("returned .. %d\n",i);
    }

    fprintf(fp_angles_ns,"\n");
    fprintf(fp_angles_ew,"\n");

#endif // ANGLE3CENTROID_M

#ifdef ANGLE3CENTROID_DIMERALPHA_M
    debug("entering ANGLE3CENTROID_DIMERALPHA_M !!!!\n");

    for (int j=0; j<dimers_alpha.size(); j++) {
        // printf("angle3dimer: %d ---> %d\n",j,dimers_alpha[j]);
        topo_sort_dimer_contacts_initially(chain_ref,dimers_alpha[j],fp_dimercontacts);
    }

#endif // ANGLE3CENTROID_DIMERALPHA_M


    /* ---------------------------------------------------------
       Microtubules (2 of 2)
       --------------------------------------------------------- */
#ifdef MTCON2
    printf("deprecated.\n");
    exit(0);


    // use chain_ref to begin chain evaluation. BEGIN NOW.
    // 1. copy int LonN, LonS, LatE, LatW to chain
    // 2.


    printf("chain_ref map of networked/connected monomers: set map for chain_later, get contacts......\n");
    for ( int i=0; i<chains_to_use; i++ ) {
        // N E S W | 1 2 3 4

        // printf(">> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);
        chain_later[i].LonN = chain_ref[i].LonN;
        chain_later[i].LatE = chain_ref[i].LatE;
        chain_later[i].LonS = chain_ref[i].LonS;
        chain_later[i].LatW = chain_ref[i].LatW;
        // printf("R> %d %d %d %d\n",chain_ref[i].LonN,chain_ref[i].LatE,chain_ref[i].LonS,chain_ref[i].LatW);
        // printf("L> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);

        // exit(0);

        Contact map1[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // print_contact_map(map1);
        build_contact_map_inter_monomer(chain_later,i,map1,1);
        // print_contact_map(map1);
        // printf("\n\n\n^^sep_vv\n\n\n");
        // print_contact_map(chain_later[i].contactsLonN);
        memcpy(chain_later[i].contactsLonN, map1, sizeof(chain_later[i].contactsLonN));
        // print_contact_map(chain_later[i].contactsLonN);
        // exit(0);

        Contact map2[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter_monomer(chain_later,i,map2,2);
        memcpy(chain_later[i].contactsLatE, map2, sizeof(chain_later[i].contactsLatE));
        // print_contact_map(chain_later[i].contactsLatE);

        Contact map3[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter_monomer(chain_later,i,map3,3);
        memcpy(chain_later[i].contactsLonS, map3, sizeof(chain_later[i].contactsLonS));
        // print_contact_map(chain_later[i].contactsLonS);

        Contact map4[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        build_contact_map_inter_monomer(chain_later,i,map4,4);
        memcpy(chain_later[i].contactsLatW, map4, sizeof(chain_later[i].contactsLatW));
        // print_contact_map(chain_later[i].contactsLatW);

        // exit(0);
    }




    // break;
    // exit(0);

    printf("> comparing contact maps: chain_ref vs. chain_later .......\n");

    for ( int i=0; i<chains_to_use; i++ ) {
    // for ( int i=183; i<184; i++ ) {
        // N E S W | 1 2 3 4
        compare_contacts(chain_ref,chain_later,i,1);
        compare_contacts(chain_ref,chain_later,i,2);
        compare_contacts(chain_ref,chain_later,i,3);
        compare_contacts(chain_ref,chain_later,i,4);
        // printf("new_totals: (NESW) %d %d %d %d\n",chain_later[i].LonN_total2x,
        //        chain_later[i].LatE_total2x,chain_later[i].LonS_total2x,chain_later[i].LatW_total2x);

    }


    // printf("201 the end.\n");
    // compare_contacts(chain_ref,chain_later,201,3);
    // compare_contacts(chain_ref,chain_later,201,1);
    // compare_contacts(chain_ref,chain_later,201,3);
    // compare_contacts(chain_ref,chain_later,201,4);

    // printf("183\n");
    // compare_contacts(chain_ref,chain_later,183,1);
    // compare_contacts(chain_ref,chain_later,183,2);
    // compare_contacts(chain_ref,chain_later,183,3);
    // compare_contacts(chain_ref,chain_later,183,4);


    // fprintf!
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,1);
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,2);
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,3);
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,4);

    // break;

#endif // END MTCON2

#ifdef MTCON_OC
    printf("deprecated.\n");
    exit(0);


    // use chain_ref to begin chain evaluation. BEGIN NOW.
    // 1. copy int LonN, LonS, LatE, LatW to chain
    // 2.

    // Contact map[MOLECULE_INITIAL_SIZE][MAX_CONTACTS];
    // build_contact_map_inter(chain_ref,1,0,map); // important: 1,0 small/big
    // memcpy(chain_ref[0].intercontacts, map, sizeof(chain_ref[0].intercontacts));



    // printf("chain_ref map of networked/connected monomers: set map for chain_later, get contacts......\n");
    // for ( int i=0; i<chains_to_use; i++ ) {
        // N E S W | 1 2 3 4
        // printf(">> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);
        // chain_later[i].LonN = chain_ref[i].LonN;
        // chain_later[i].LatE = chain_ref[i].LatE;
        // chain_later[i].LonS = chain_ref[i].LonS;
        // chain_later[i].LatW = chain_ref[i].LatW;
        // printf("R> %d %d %d %d\n",chain_ref[i].LonN,chain_ref[i].LatE,chain_ref[i].LonS,chain_ref[i].LatW);
        // printf("L> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);

        // exit(0);

        // Contact map1[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // print_contact_map(map1);
        // build_contact_map_inter_monomer(chain_later,i,map1,1);
        // print_contact_map(map1);
        // printf("\n\n\n^^sep_vv\n\n\n");
        // print_contact_map(chain_later[i].contactsLonN);
        // memcpy(chain_later[i].contactsLonN, map1, sizeof(chain_later[i].contactsLonN));
        // print_contact_map(chain_later[i].contactsLonN);
        // exit(0);

        // Contact map2[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // build_contact_map_inter_monomer(chain_later,i,map2,2);
        // memcpy(chain_later[i].contactsLatE, map2, sizeof(chain_later[i].contactsLatE));
        // print_contact_map(chain_later[i].contactsLatE);

        // Contact map3[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // build_contact_map_inter_monomer(chain_later,i,map3,3);
        // memcpy(chain_later[i].contactsLonS, map3, sizeof(chain_later[i].contactsLonS));
        // print_contact_map(chain_later[i].contactsLonS);

        // Contact map4[MOLECULE_INITIAL_SIZE_MT][MAX_CONTACTS];
        // build_contact_map_inter_monomer(chain_later,i,map4,4);
        // memcpy(chain_later[i].contactsLatW, map4, sizeof(chain_later[i].contactsLatW));
        // print_contact_map(chain_later[i].contactsLatW);

        // exit(0);
    // }

    printf("> comparing contact maps: chain_ref vs. chain_later .......\n");
    for ( int i=0; i<chains_to_use; i++ ) {
    // for ( int i=183; i<184; i++ ) {
        // N E S W | 1 2 3 4

        evaluate_original_contacts_now(chain_ref,chain_later,i,1);
        evaluate_original_contacts_now(chain_ref,chain_later,i,2);
        evaluate_original_contacts_now(chain_ref,chain_later,i,3);
        evaluate_original_contacts_now(chain_ref,chain_later,i,4);
        // compare_contacts(chain_ref,chain_later,i,1);
        // compare_contacts(chain_ref,chain_later,i,2);
        // compare_contacts(chain_ref,chain_later,i,3);
        // compare_contacts(chain_ref,chain_later,i,4);
        // printf("new_totals: (NESW) %d %d %d %d\n",chain_later[i].LonN_total2x,
        //        chain_later[i].LatE_total2x,chain_later[i].LonS_total2x,chain_later[i].LatW_total2x);
    }

    // fprintf!
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,1);
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,2);
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,3);
    fprintf_contacts_for_one_monomer_face(chain_later,chains_to_use,4);
#endif // END MTCON_OC


#ifdef DIMERMAP_2
    printf("evaluating section dimermap_2 now ..\n");

    // angles:
    for (int i=0; i<dimermap_alphas.size(); i++)
    {
        maprow_angles = get_3angles_dimerbyalpha(chain_ref,chain_later,dimermap_alphas[i]);
        map_angles.insert( std::pair<int, std::vector<double> >(i,maprow_angles));
        maprow_angles.clear();
    }


    // for ( int i=0; i < chains_to_use; i++)
    // {
    //     chain_later[i].ComputeCentroid();
    // }

    for (int i=0; i<chains_to_use; i++ )
    {
        // N E S W | 1 2 3 4
        debug(">> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);
        chain_later[i].num_atoms_ca = chain_ref[i].num_atoms_ca;
        chain_later[i].LonN = chain_ref[i].LonN;
        chain_later[i].LatE = chain_ref[i].LatE;
        chain_later[i].LonS = chain_ref[i].LonS;
        chain_later[i].LatW = chain_ref[i].LatW;
        debug("R> %d %d %d %d\n",chain_ref[i].LonN,chain_ref[i].LatE,chain_ref[i].LonS,chain_ref[i].LatW);
        debug("L> %d %d %d %d\n",chain_later[i].LonN,chain_later[i].LatE,chain_later[i].LonS,chain_later[i].LatW);

        debug("pos_later-35: %f %f %f\n",chain_later[i].pos[34].x, \
              chain_later[i].pos[34].y, \
              chain_later[i].pos[34].z);
    }

    // ComputeCentroid()
    for (int i=0; i < chains_to_use; i++)
    {
        // chain_later[i].print_centroid();
        chain_later[i].ComputeCentroid();
        // chain_later[i].print_centroid();
    }

    // // curvature:
    // Vector vec_zero;
    // vec_zero.x = vec_zero.y = vec_zero.z = 0.0;

    for (int i=0; i<dimermap_alphas.size(); i++)
    {
        int a_j, b_j;
        Curvature curv_ai, curv_bi, curv_long;

        // alat
        a_j = dimermap_alphas[i];
        b_j = chain_ref[a_j].LonN;
        // b_j = a_j + 1;

        int a_west,b_west;
        int a_east,b_east;
        int b_north;
        a_west = b_west = a_east = b_east = b_north = 0;

        a_west = chain_ref[a_j].LatW;
        a_east = chain_ref[a_j].LatE;
        b_west = chain_ref[b_j].LatW;
        b_east = chain_ref[b_j].LatE;
        b_north= chain_ref[b_j].LonN;

        // printf("a-b: %d %d\n",a_j,b_j);
        // printf("awe,bwen: %d %d   %d %d %d\n",a_west,a_east,\
        //        b_west,b_east,b_north);
        // chain_later[a_j].print_centroid();
        // chain_later[b_j].print_centroid();

        // FAIL
        // printf("a_j centroid: %f %f %f\n",chain_later[a_j].centroid); // FAIL, cen.x,.y,.z
        // printf("b_j centroid: %f %f %f\n",chain_later[b_j].centroid);




        //     (chain_ref[a_j].LatW == -1) ||
        //     (chain_ref[a_j].LatE == -1) ||
        //     (chain_ref[b_j].LatW == -1) ||
        //     (chain_ref[b_j].LatE == -1) ||
        //     (chain_ref[b_j].LonN == -1))

        // if (
        //     (chain_ref[a_j].LatW == -1) ||
        //     (chain_ref[a_j].LatE == -1) ||
        //     (chain_ref[b_j].LatW == -1) ||
        //     (chain_ref[b_j].LatE == -1) ||
        //     (chain_ref[b_j].LonN == -1))
        // {
        //     printf("alpha: %d\n",a_j);
        //     printf("l %d   %d r\n",chain_ref[a_j].LatW,chain_ref[a_j].LatE);
        //     maprow_curv.push_back(0.0);
        //     continue;
        // }
        // else
        // {


        // TESTING
        // if((a_j != -1) && (a_west != -1) && (a_east != -1))
        // {
        //     printf("NO-MINUS-ONEs  awe:  %d %d %d\n",a_j,a_west,a_east);
        // }
        // else
        // {
        //     printf("MINUS          awe:  %d %d %d\n",a_j,a_west,a_east);
        // }


        if((a_j != -1) && (a_west != -1) && (a_east != -1))
        {
            // printf("N0-MINUS-ONE  a-b: %d %d\n",a_j,b_j);
            // printf("awe,bwen: %d %d   %d %d %d\n",a_west,a_east,\
            //        b_west,b_east,b_north);

            curv_ai.p1 = chain_later[a_west].centroid;
            curv_ai.p2 = chain_later[a_j].centroid;
            curv_ai.p3 = chain_later[a_east].centroid;
            // curv_ai.print_points();
            curv_ai.get_radius_circumscribed();
        }
        else
        {
            curv_ai.radius = 0.0;
        }
        // curv_ai.radius = 0.0;
        maprow_curv.push_back(curv_ai.radius);


        if((b_j != -1) && (b_west != -1) && (b_east != -1))
        {
            curv_bi.p1 = chain_later[b_west].centroid;
            curv_bi.p2 = chain_later[b_j].centroid;
            curv_bi.p3 = chain_later[b_east].centroid;
            // curv_bi.print_points();
            curv_bi.get_radius_circumscribed();
        }
        else
        {
            curv_bi.radius = 0.0;
        }
        // curv_bi.radius = 0.0;
        maprow_curv.push_back(curv_bi.radius);



        if (b_north != -1)
        {
            curv_long.p1 = chain_later[a_j].centroid;
            curv_long.p2 = chain_later[b_j].centroid;
            curv_long.p3 = chain_later[b_north].centroid;
            // curv_long.print_points();
            curv_long.get_radius_circumscribed();
        }
        else
        {
            curv_long.radius = 0.0;
        }
        // curv_long.radius = 0.0;
        maprow_curv.push_back(curv_long.radius);


        // curv_ai.p1 = chain_ref[a_j].centroid;
        // curv_ai.p2 = chain_ref[chain_ref[a_j].LatW].centroid;
        // curv_ai.p3 = chain_ref[chain_ref[a_j].LatE].centroid;
        // // curv_ai.print_points();
        // // printf("\t-->\n");
        // // curv_ai.points_ascending_z();
        // debug("\ncurv_ai: %d\n",a_j);
        // curv_ai.print_points();
        // // curv_ai.get_circumscribed_circle();
        // curv_ai.get_radius_circumscribed();
        // debug(">: %f %f\n",curv_ai.radius,curv_ai.curvature);
        // // printf("alat ^\n");
        // maprow_curv.push_back(curv_ai.radius);


        // // blat
        // curv_bi.p1 = chain_ref[b_j].centroid;
        // curv_bi.p2 = chain_ref[chain_ref[b_j].LatW].centroid;
        // curv_bi.p3 = chain_ref[chain_ref[b_j].LatE].centroid;
        // // curv_bi.print_points();
        // // printf("\t-->\n");
        // // curv_bi.points_ascending_z();
        // // curv_bi.print_points();
        // // curv_bi.get_circumscribed_circle();
        // curv_bi.get_radius_circumscribed();
        // // printf(">: %f %f\n",curv_bi.radius,curv_bi.curvature);
        // // printf("blat ^\n");
        // maprow_curv.push_back(curv_bi.radius);



        // // blon
        // curv_long.p1 = chain_ref[a_j].centroid;
        // curv_long.p2 = chain_ref[b_j].centroid;
        // curv_long.p3 = chain_ref[chain_ref[b_j].LonN].centroid;
        // // curv_long.print_points();
        // // printf("\t-->\n");
        // // curv_long.points_ascending_z();
        // // curv_long.print_points();
        // // curv_long.get_circumscribed_circle();
        // curv_long.get_radius_circumscribed();
        // // printf(">: %f %f\n",curv_long.radius,curv_long.curvature);
        // // printf("blon ^\n");
        // maprow_curv.push_back(curv_long.radius);
        // // }


        map_curv.insert(std::pair<int, std::vector<double> >(i,maprow_curv));
        // map_contacts.insert( std::pair<int, std::vector<int> >(i,maprow_contacts));

        maprow_curv.clear();

    }
    // exit(0


    // printf("contact_check: %d\n",map_contacts[6][3]);


#ifndef NDEBUG
    // print map, check.
    printf("contacts1:\n");
    for (int j=0; j<map_contacts.size(); j++)
    {
        printf("\t%d: %d",j,map_contacts[j].size());
        for (int k=0; k<map_contacts[j].size(); k++)
        {
            printf("  %d-%d",k,map_contacts[j][k]);
        }
        printf("\n");
    }
#endif


#ifndef NDEBUG
    // print map, check.
    printf("angles:\n");
    for (int j=0; j<map_angles.size(); j++)
    {
        printf("\t%d: %d",j,map_angles[j].size());
        for (int k=0; k<map_angles[j].size(); k++)
        {
            printf("   %d %f",k,map_angles[j][k]);
        }
        printf("\n");
    }
#endif

    // printf("angles:\n");
    // for (int j=0; j<map_angles.size(); j++)
    // {
    //     printf("\t%d:%d->",j,map_angles[j]);
    //     for (int k=0; k<map_angles[j].size(); k++)
    //     {
    //         printf(" %d",j,map_angles[j][k]);
    //     }
    //     printf("\n");
    // }

    // printf("curv:\n");
    // for (int j=0; j<map_curv.size(); j++)
    // {
    //     printf("\t%d:%d->",j,map_curv[j]);
    //     for (int k=0; k<map_curv[j].size(); k++)
    //     {
    //         printf(" %d",j,map_curv[j][k]);
    //     }
    //     printf("\n");
    // }



    // max:
    fprintf(fp_mt_analysis,"# frame: %d\n",frame_position);

    // write data!
    dimer_east = dimer_west = dimer_ext = alpha = beta = 0;

    for (int i=0; i<dimermap_alphas.size(); i++)
    {
        alpha = dimermap_alphas[i];
        beta = chain_ref[alpha].LonN;
        dimer_east = map_contacts[alpha][1] + map_contacts[beta][1];
        dimer_west = map_contacts[alpha][3] + map_contacts[beta][3];
        dimer_ext = dimer_east + dimer_west + map_contacts[alpha][2] + map_contacts[beta][0];


        // S,I,N,E,W,E
        fprintf(fp_mt_analysis,"%3d %6d %6d %3d %6d %6d  %3d %3d %3d %3d %3d %3d", \
                alpha,chain_ref[alpha].index,chain_ref[alpha].findex,   \
                beta,chain_ref[beta].index,chain_ref[beta].findex,      \
                map_contacts[alpha][2],map_contacts[alpha][0],          \
                map_contacts[beta][0],dimer_east,dimer_west,dimer_ext);

        // alpha lateral, beta lateral, beta longitudinal
        fprintf(fp_mt_analysis," %5.1f %5.1f %5.1f",                  \
                map_angles[i][0],map_angles[i][1],map_angles[i][2]);


        // alat, blat, blon
        fprintf(fp_mt_analysis,"  %7.1f %7.1f %7.1f\n", \
                map_curv[i][0], \
                map_curv[i][1], \
                map_curv[i][2]);
    }
    // fprintf(fp_mt_analysis,"\n");



    // map clear.
    map_contacts.clear();
    map_angles.clear();
    map_curv.clear();

#endif // DIMERMAP_2


#ifdef DCD_WRITE
    // READ
    // static void *open_dcd_read(const char *path, const char *filetype,
    //                            int *natoms) {
    // int rc = read_next_timestep(v, natoms, timestep);

    // WRITE
    // static void *open_dcd_write(const char *path, const char *filetype,
    //                             int natoms) {
    // static void close_file_write(void *v) {

    // open_dcd_write(fn_dcd_write,"dcd",&natoms);

    // write_dcdstep
    // static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N,
    //                          const float *X, const float *Y, const float *Z,
    //                          const double *unitcell, int charmm) {

    // static int write_dcdheader(fio_fd fd, const char *remarks, int N,
    //                            int ISTART, int NSAVC, double DELTA, int with_unitcell,
    //                            int charmm) {

    // load_chain_coords_to_timestep(dcdw,chain_later,chains_to_use,vw,&timestep_w,natoms_w);
    // load_chain_coords_to_timestep(chain_later,num_chains,&timestep_w);
    load_chain_to_timestep(chain_later,num_chains,&timestep_w);


#ifdef DCD_WRITE_UNMOD
    // Write the DCD read in.
    write_timestep(vw,&timestep);
// #elif DCD_WRITE_ROT

#else
    // Write modified coordinates.
    write_timestep(vw,&timestep_w);
#endif // A
#endif


    // delete
#ifdef DEFAULT
    delete [] chain_ref;
    delete [] chain_later;
#endif
#ifdef ONEMOL
    delete [] chain_ref;
#endif

#ifdef DCDREAD
    // } // DCD PRIMARY LOOP

        debug("current: %d\n",nset2);
        printf("frame: --> %d <-- was evaluated.\n",frame_position);


#ifndef NDEBUG
        printf("ref-findex(%d): %f\n",chain_ref[0].findex,chain_ref[0].pos[chain_ref[0].findex].y);
        printf("0-findex(%d): %f\n",chain_0[0].findex,chain_0[0].pos[chain_0[0].findex].y);
        printf("later-findex(%d): %f\n",chain_later[0].findex,chain_later[0].pos[chain_later[0].findex].y);
        if (chains_to_use > 1) {
            printf("later-findex(%d): %f\n",chain_later[1].findex,chain_later[1].pos[chain_later[1].findex-chain_later[0].num_atoms_ca].y);
        }
#endif


        if (nset2 + advance_size + 1 <= stop) {
            frame_position += advance_dcd(dcd->nsets,advance_size,dcd,natoms,&timestep);
            printf("frame: --> %d <-- loaded.\n",frame_position);
            load_dcd_to_chain(dcd,chain_later,num_chains);
            // nset2 += advance_size + 1;
        }
        nset2 += advance_size + 1;
    } while (nset2<=stop);

    debug("..closing dcd..\n");
    close_file_read(v);
    // END DCD
    // typedef struct {
    //     fio_fd fd;
    //     int natoms; 382
    //     int nsets; 17501 (0-17500)
    //     int setsread;
    //     int istart;
    //     int nsavc;
    //     double delta;
    //     int nfixed;
    //     float *x, *y, *z; ->x[0-381];
    //     int *freeind;
    //     float *fixedcoords;
    //     int reverse;
    //     int charmm;
    //     int first;
    //     int with_unitcell;
    // } dcdhandle;

    delete [] chain_ref;
    // delete [] chain;
    delete [] chain_0;
    delete [] chain_later;

    printf("DCDREAD complete.\n\tThe maximum possible frame_position was: %d\n",stop);
    printf("\tThe last frame evaluated was: %d\n",frame_position);
#endif //DCDREAD

#ifdef DCD_WRITE_E

    // open_dcd_write(fn_dcd_write,"dcd",natoms);
    // static void close_file_write(void *v) {
    close_file_write(vw);

#endif // DCD_WRITE_E


#ifdef DIMERMAP_3
    printf("dimermap_3 completing now ..\n");

    fclose(fp_mt_analysis);

#endif // DIMERMAP_3

#ifdef CURVATURE_E
    // CURVATURE_E
    // FILE * fp_proto_centroid_dist;
    // fp_proto_centroid_dist = fopen("distance_proto_centroids.dat", "a+");


    // last evaluated.
    fprintf(fp_proto_centroid_dist,"# last_frame: %d\n",frame_position);
    fprintf(fp_curvature_global,"# last_frame: %d\n",frame_position);
    fprintf(fp_radius_of_curvature_global,"# last_frame: %d\n",frame_position);
    fprintf(fp_curvature_local,"# last_frame: %d\n",frame_position);
    fprintf(fp_radius_of_curvature_local,"# last_frame: %d\n",frame_position);


    // stop
    fprintf(fp_proto_centroid_dist,"# stop_frame: %d\n",stop);
    fprintf(fp_curvature_global,"# stop_frame: %d\n",stop);
    fprintf(fp_radius_of_curvature_global,"# stop_frame: %d\n",stop);
    fprintf(fp_curvature_local,"# stop_frame: %d\n",stop);
    fprintf(fp_radius_of_curvature_local,"# stop_frame: %d\n",stop);

    // centroid distance.
    fclose(fp_proto_centroid_dist);
    // global.
    fclose(fp_curvature_global);
    fclose(fp_radius_of_curvature_global);
    // local.
    fclose(fp_curvature_local);
    fclose(fp_radius_of_curvature_local);

#endif // CURVATURE_E

#ifdef CHI_END

    fclose(fp_chi);

#endif // CHI_END

#ifdef ANGLE3CENTROID_E
    // fprintf(fp_angles_ns, "# time later Center|West|East or C|South|North* (frame:%d)\n",frame_position);
    // fprintf(fp_angles_ew, "# time later Center|West|East* or C|South|North (frame:%d)\n",frame_position);


    // last
    fprintf(fp_angles_ns,"# last_frame: %d\n",frame_position);
    fprintf(fp_angles_ew,"# last_frame: %d\n",frame_position);
    // stop
    fprintf(fp_angles_ns,"# stop_frame: %d\n",stop);
    fprintf(fp_angles_ew,"# stop_frame: %d\n",stop);


    printf("Angle analysis by the position of 3 centroids is complete.\n");
    fclose(fp_angles_ns);
    fclose(fp_angles_ew);
#endif // ANGLE3CENTROID_E

#ifdef CENTROIDMOVEMENT_E

    printf("Centroid movement by position and net distance tracking is complete.\n");
    // stop
    fprintf(fp_xyzpos,"# last_frame: %d\n",frame_position);
    fprintf(fp_xyzpos,"# stop_frame: %d\n",stop);
    fclose(fp_xyzpos);

#endif // CENTROIDMOVEMENT_E


#ifdef CONTACTMAP_SM3
    // FILE * fp_contacts
    // fp_contacts_n = fopen("contacts_intra.dat", "w+");
    fclose(fp_contacts);


#ifdef CONTACTSBYRESFILE1
    fclose(fp_contactsresfile1);
#endif // CONTACTSBYRESFILE1

#endif // CONTACTMAP_SM3

#ifdef CONTACTMAP_MT3
    // if(chains_to_use > 1){

#ifdef CONFILE
    fclose(fp_contacts_n);
    fclose(fp_contacts_e);
    fclose(fp_contacts_s);
    fclose(fp_contacts_w);
#endif // CONFILE

    fclose(fp_dimercontacts);
    // }
#endif // CONTACTMAP_MT3

#ifdef MTPF_E
printf("MTPF_E: processing..\n");

for(int i = 0; i < pf_array.size()*4; i++){
    if (outcontactfile[i].is_open()) {
        outcontactfile[i].close();
    }
 }


printf("MTPF_E: complete.\n");
#endif // MTPF_E


    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    std::cout << "\nexiting successfully ...\n";
    return 0;
}
