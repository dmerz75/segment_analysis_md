// readfile.cpp

/* ---------------------------------------------------------
   standard libraries
   --------------------------------------------------------- */
// #include <stdio.h>
#include <stdlib.h> // strtod?, stod
// #include <assert.h>
#include <string.h>
// #include <string>
#include <iostream>
#include <string>


#include <fstream>
#include <stdlib.h>
#include <stdio.h>


#include <algorithm> // count
#include <iterator> // istream_iterator
#include <ctime>

// namespaces
// using namespace

// extern
extern int endoffile;

/* ---------------------------------------------------------
   headers
   --------------------------------------------------------- */
#include "readfile.h"

#include "chain.h"
// #define MOLECULE_INITIAL_SIZE 700
// Chain chain;



// Finding symbol: ReadLines
// Database directory: /home/dmerz3/sop_dev/contacts/segment/

// *** include/readfile.h:
// <global>[27]                   void ReadLines (Chain *chain_segment);

// *** src/main.cpp:

// From DCDREAD
// main[102]                      ReadLines(&chain_pdb[i]);
// main[153]                      ReadLines(&chain_ref[i]);
// main[194]                      ReadLines(&chain_0[i]);
// main[236]                      ReadLines(&chain_later[i]);

// From DEFAULT
// main[1189]                     ReadLines(&chain_ref[i]);
// main[1225]                     ReadLines(&chain_later[i]);

// From ONEMOL
// main[1285]                     ReadLines(&chain_ref[i]);

// *** src/readfile.cpp:
// ReadLines[43]                  void ReadLines (Chain *chain) {

//     Search complete.  Search time = 0.07 seconds.



/* ---------------------------------------------------------
   functions
   --------------------------------------------------------- */
void ReadLines (Chain *chain) {
    double x, y, z;
    int count_line=0, count_ca=0, count_index=0;

    // PRINT HERE
    // std::cout << "readcheck:" << std::endl;
    // std::cout << "line: " << count_line << " CA: " << count_ca << " index: " << count_index << std::endl;
    // std::cout << "begin: " << chain->file_line_begin << " end: " << chain->file_line_end << std::endl;
    // std::cout << "index: " << chain->index << " findex: " << chain->findex << std::endl;


    // ctime (1 of 2)
    clock_t begin = clock();

    // get file, open it
    // std::fstream readthisfile (fname);
    std::fstream readthisfile (chain->filename);
    std::string line;
    // std::string pos_unit;
    // char pdb_line[83];

    if (readthisfile.is_open()) {
        // call read/search function
#ifdef INFO
        std::cout << "reading PDB ...\n";

#endif

        // loop
        while (std::getline(readthisfile,line)) {

            count_line ++;
            if ( count_line <= chain->file_line_begin ) {
                // std::cout << count_line << " vs. " << chain->file_line_begin << std::endl;
                // std::cout << line << '\n';
                continue;
            }
            // PRINT HERE
            // std::cout << line << '\n';

            // printf("%s\n",line);

            // PRINT HERE
            // std::cout << line << '\n';
            std::string begin4 = line.substr(0,4);
            std::size_t foundCA = line.find("CA");

            // usually prints "13" because CA is there.
            // if (foundCA!=std::string::npos) {
            //     std::cout << foundCA << std::endl;
            // }

            // if the first 4 positions == ATOM  and
            // line.find("CA") finds something (see next line, just below)


            if ((begin4 == "ATOM") && (foundCA!=std::string::npos)) {

                count_ca ++;
                count_index ++;
                int i;
                i = count_index - 1;

                // std::cout << "atom_pos: " << '\t' << std::endl;
                std::string str_x = line.substr(30,8);
                std::string str_y = line.substr(38,8);
                std::string str_z = line.substr(46,8);
                // std::cout << str_x << str_y << str_z << '\n' << std::endl;

                // std::cout << "x: " << chain->pos(i,0) << " y: " << chain->pos(i,1) << " z: " << chain->pos(i,2) << std::endl;
                // amazing!! str_x.c_str()
                x = atof(str_x.c_str());
                y = atof(str_y.c_str());
                z = atof(str_z.c_str());

                // armadillo (arma) namespace
                // chain->pos(i,0) = x;
                // chain->pos(i,1) = y;
                // chain->pos(i,2) = z;
                // std::cout <<
                // std::cout << "x: " << chain->pos(i,0) << " y: " << chain->pos(i,1) << " z: " << chain->pos(i,2) << std::endl;

                // Vector pos[MOLECULE_INIT_SIZE];
                chain->pos[i].x = x;
                chain->pos[i].y = y;
                chain->pos[i].z = z;

                // PRINT HERE
                // std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;

                // get resname, resid, indices;
                // memcpy(resname, &str_dest[17],3);
                // resname[3] = '\0';
                // strcpy(mol->resname[count_ATOM_CA],resname);
                // memcpy(resid, &str_dest[23],3);


                // currently, resname is blank.
                // std::cout << "the value of the resname is: " << chain->resname[count_ca-1] << std::endl;
                std::string str_resname = line.substr(17,3);
                strcpy(chain->resname[count_ca-1],str_resname.c_str());
                // chain->resname[4] = '\0';
                // std::cout << "the value of the resname is: " << chain->resname[count_ca-1] << std::endl;
                chain->resid[count_ca-1] = atoi(line.substr(23,3).c_str());
                // std::cout << "the value of the resid is: " << chain->resid[count_ca-1] << std::endl;


                // first breakout.
                if (count_ca > 1) {
                    // PRINT HERE
                    // std::cout << chain->resid[count_ca-2] << " vs. " << chain->resid[count_ca-1] << std::endl;
                    if (chain->resid[count_ca-2] >= chain->resid[count_ca-1]) {


                        // chain->file_line_end = chain->file_line_begin + count_line - 1;
                        chain->file_line_end = count_line - 1;
                        // chain->findex = chain->index + count_ca - 1;


                        // PRINT HERE
                        // std::cout << "readcheck:" << std::endl;
                        // std::cout << "line: " << count_line << " CA: " << count_ca - 1 << " index: " << count_index << std::endl;
                        // std::cout << "begin: " << chain->file_line_begin << " end: " << chain->file_line_end << std::endl;
                        // std::cout << "index: " << chain->index << " findex: " << chain->findex << std::endl;



#ifdef INFO
                        // ctime (2 of 2)
                        clock_t end = clock();
                        double elapsed_time = double(end - begin) / CLOCKS_PER_SEC; // only second accuracy
                        printf("read some of file in %f seconds.\n",elapsed_time);
#endif //
                        chain->num_atoms_ca = count_ca - 1;
                        return;
                    }
                }
            } // end of ATOM && CA
        } // while getline
        // End of file.

        // chain->file_line_end = chain->file_line_begin + count_line;
        chain->file_line_end = count_line - 1;
        chain->num_atoms_ca = count_ca;
        // chain->findex = chain->index + count_ca -1;

        // ctime (2 of 2)
        clock_t end = clock();
        double elapsed_time = double(end - begin) / CLOCKS_PER_SEC; // only second accuracy

#ifdef INFO
        printf("read to the end of the file in %f seconds.\n",elapsed_time);
        std::cout << "THE FILE LINE COUNTS begin==end? " << chain->file_line_begin << " ?= " << chain->file_line_end << std::endl;
#endif

        std::cout << "\tnum_chains_arg should be " << chain->chainid + 1 << "\n"<<std::endl;
    } else {
        printf("Error opening file!\n");
        printf("could not open %s\n",chain->filename);
        // std::fstream readthisfile (chain->filename);
        exit (EXIT_FAILURE);
    }
    return;
}

void ReadEveryLine (Chain *chain) {
    double x, y, z;
    int count_line=0, count_ca= 0, count_atom=0, count_index=0, count_ca_index;
    int i = 0, whilecheck=0;

    // ctime (1 of 2)
    clock_t begin = clock();

    std::fstream readthisfile (chain->filename);
    std::string line;

    if (readthisfile.is_open()) {

#ifdef INFO
        std::cout << "\nreading PDB ...  [starting line:(" \
                  << chain->file_line_begin << ")]. \n" << std::endl;
#endif

        // loop
        while (std::getline(readthisfile,line)) {

            count_line ++;
            if ( count_line <= chain->file_line_begin ) {
                continue;
            }
            // // lines actually read during each iteration.
            // printf("%s\n",line.c_str());

            std::string begin3 = line.substr(0,3);
            std::string begin4 = line.substr(0,4);
            std::size_t foundCA = line.find("CA");

            // ATOM, HETATM
            if ((begin4 == "ATOM") || (begin4 == "HETA")) {

                // if ((begin4 == "ATOM") && (foundCA!=std::string::npos)) {
                if (foundCA!=std::string::npos) {
                    count_ca++;
                    count_ca_index++;
                }

                count_atom ++;
                count_index ++;
                i = count_index - 1;

                // ATOM   4698  CB  PRO A 558       6.662 -10.694 -99.521  1.00 99.99           C
                std::string str_atomtype = line.substr(12,3);
                std::string str_resname = line.substr(17,3);
                std::string str_resid = line.substr(23,3);
                std::string str_x = line.substr(30,8);
                std::string str_y = line.substr(38,8);
                std::string str_z = line.substr(46,8);

                // resname,resid,x,y,z;
                strcpy(chain->atomtype[i],str_atomtype.c_str());
                strcpy(chain->resname[i],str_resname.c_str());
                chain->resid[i] = atoi(str_resid.c_str());
                x = atof(str_x.c_str());
                y = atof(str_y.c_str());
                z = atof(str_z.c_str());
                chain->pos[i].x = x;
                chain->pos[i].y = y;
                chain->pos[i].z = z;

#ifdef INFO
                // ctime (2 of 2)
                clock_t end = clock();
                double elapsed_time = double(end - begin) / CLOCKS_PER_SEC; // only second accuracy
                printf("read some of file in %f seconds.\n",elapsed_time);


                // For any "ATOM","HETA"
                printf("previous: %d   current: %d   Atoms: %d     ca: %d \n", \
                       chain->resid[i-1],chain->resid[i],               \
                       count_atom,count_ca);

                // // For any "ATOM","HETA"
                // printf("previous: %d   current: %d   Atoms: %d     ca: %d \n", \
                //        chain->resid[i-1],chain->resid[i],               \
                //        chain->num_atoms,chain->num_atoms_ca);
#endif //

                if ((chain->resid[i-1] != chain->resid[i]) &&         \
                    (chain->resid[i-1] + 1 != chain->resid[i]) && \
                    (chain->resid[i-1] != 0)) {

                    // need the minus 1.
                    chain->file_line_end = count_line - 1;

                    // if ((begin4 == "ATOM") || (begin4 == "HETA")) {
                    chain->num_atoms = count_atom - 1;

                    if (foundCA!=std::string::npos) {
                        chain->num_atoms_ca = count_ca - 1;
                    } else {
                        chain->num_atoms_ca = count_ca;
                    }

                    // printf("RESID BREAK FOUND!\n");

                    // ctime (2 of 2)
                    clock_t end = clock();
                    double elapsed_time = double(end - begin) / CLOCKS_PER_SEC; // only second accuracy
#ifdef INFO
                    printf("read to the end of the file in %f seconds.\n",elapsed_time);
                    std::cout << "THE FILE LINE COUNTS begin==end:" << chain->file_line_begin << " <=> " << chain->file_line_end << std::endl;
#endif
                    std::cout << "\tnum_chains_arg should be " << chain->chainid + 1 << "\n"<<std::endl;
                    return;
                } // if break chain!!!
            } // end of ATOM in first4
        } // while getline

        // End of file.
        chain->file_line_end = count_line;
        chain->num_atoms_ca = count_ca;
        chain->num_atoms = count_atom;

        // printf("ATOMS: %d  CALPHAs: %d\n",chain->num_atoms,chain->num_atoms_ca);
        // ctime (2 of 2)

        clock_t end = clock();
        double elapsed_time = double(end - begin) / CLOCKS_PER_SEC; // only second accuracy
// #ifdef INFO
        printf("read to the end of the file in %f seconds.\n",elapsed_time);
        std::cout << "THE FILE LINE COUNTS begin==end:" << chain->file_line_begin << " <=> " << chain->file_line_end << std::endl;
// #endif
        std::cout << "\tnum_chains_arg should be " << chain->chainid + 1 << "\n"<<std::endl;

        return;

    } else { // if file did not open.
        printf("Error opening file!\n");
        printf("could not open %s\n",chain->filename);
        // std::fstream readthisfile (chain->filename);
        exit (EXIT_FAILURE);
    } // if file is open
    return;
}

void ReadMolecularContent (char *fname) {
// int ReadMolecularContent ( std::string fname ) {

    // std::string pdbfilename;
    // pdbfilename = PdbFile.filename();

    // ctime (1 of 2)
    clock_t begin = clock();

    // get file, open it
    std::ifstream readthisfile (fname);
    std::string line;
    // char pdb_line[83];

    if (readthisfile.is_open()) {
        // call read/search function
        // std::cout << "reading PDB ...\n";

        // new lines will be skipped unless we stop it from happening:
        readthisfile.unsetf(std::ios_base::skipws);

        // 180135 :: 0.29-0.30s
        // count the newlines with an algorithm specialized for counting:
        unsigned line_count = count(std::istream_iterator<char>(readthisfile),
                                    std::istream_iterator<char>(),'\n');

        // 180135 :: 0.47s
        // unsigned line_count = count(istreambuf_iterator<char>(readthisfile),
        //                             istreambuf_iterator<char>(),'\n');

        // loop
        while (std::getline(readthisfile,line)) {
            // printf("%s\n",line);
            std::cout << line << '\n';
        }
        // ctime (2 of 2)
        clock_t end = clock();
        double elapsed_time = double(end - begin) / CLOCKS_PER_SEC; // only second accuracy
        printf("read %d lines in %f seconds.\n",line_count,elapsed_time);


        // readthisfile.close()
    } else {
        printf("Error opening file!\n");
        exit (EXIT_FAILURE);
    }
    return;
}
