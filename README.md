# Segment Analysis for Coarse-grained MD in the SOP model
The self-organized polymer (SOP) model is a coarse-grained molecular dynamics software that
employs Brownian motion for the dynamics describing macromolecular structures. It uses a finite-extensible
non-linear elastic potential for the covalent bonds, a Lennard-Jones potential for non-covalent native contacts,
and the repulsive part of the Lennard-Jones potential for non-covalent non-native contacts. Due to the size of the
macromolecular structures and the corresponding trajectories acquired using the SOP model, certain scripting languages
are not fast enough or capable of holding in memory the necessary information to perform routine analyses.
Therefore, Segment Analysis was written to facilitate the analysis of these structures and their dynamics.

## Quickstart:
    make all
Compiling all of the softwares available is a quick way to get started. See the [Makefile](blob/master/Makefile).
The following softwares will be compiled:
* dimermap_release
* indices_release
* topology_release
* write_dcd_release
* costension_release
* curva_mod_release
* curva_mod3_release
* anglec_release
* cenmov_release
* contactmapM_release
* contactmapP_release
* contactmapN_release
* dcd_release
* chi_release


## bin: executables
## build: all object files, 'make clean' to remove
## doc: notes, etc
## include: header files
## lib: libraries ..
## src: The application and only the application’s source files.
## test: testing

# Available Analyses from the Makefile:
contact_release
proto_release
mtcon_release
