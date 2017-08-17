# Segment Analysis for Coarse-grained MD in the SOP model
The self-organized polymer (SOP) model is a coarse-grained molecular dynamics software that
employs Brownian motion for the dynamics describing macromolecular structures. It uses a finite-extensible
non-linear elastic potential for the covalent bonds, a Lennard-Jones potential for non-covalent native contacts,
and the repulsive part of the Lennard-Jones potential for non-covalent non-native contacts. Due to the size of the
macromolecular structures and the corresponding trajectories acquired using the SOP model, certain scripting languages
are not fast enough or capable of holding in memory the necessary information to perform routine analyses.
Therefore, Segment Analysis was written to facilitate the analysis of these structures and their dynamics.

## Quickstart:
Create a bin directory if it does not yet exist. Then run:

    mkdir -p bin
    make all

After creating all of the executables, consider copying them to /usr/local/bin/.

    sudo cp bin/run_segment* /usr/local/bin

Compiling all of the softwares available is a quick way to get started. See the [Makefile](./Makefile).
The following softwares will be compiled (target:executable). $(EXEC) is usually run_segment.
* dimermap_release: bin/$(EXEC)_dcd_dimermap_mt
* indices_release:
* topology_release: bin/$(EXEC)_top
* write_dcd_release: bin/$(EXEC)_dcd_write
* costension_release:
* curva_mod_release:
* curva_mod3_release: bin/$(EXEC)_dcd_curva3_pfm
* anglec_release:
    * bin/$(EXEC)_dcd_anglec_mt
    * bin/$(EXEC)_dcd_anglec_pf
* cenmov_release: bin/$(EXEC)_dcd_cenmov_mt
* contactmapM_release: bin/$(EXEC)_dcd_contactmap_mt
* contactmapP_release: bin/$(EXEC)_dcd_contactmap_pf
* contactmapN_release: bin/$(EXEC)_dcd_contactmap_1
* dcd_release:  bin/$(EXEC)_dcd
* chi_release: bin/$(EXEC)_dcd_chi_1
