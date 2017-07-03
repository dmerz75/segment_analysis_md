# Makefile

# Declaration of variables
CC       := gcc
CXX      := g++
CL       := clang --analyze
# LINKER := gcc -fPIC
# OBJDIR := obj

# Compiler Flags: Use $(CF) for generic/old architectures
CF       := -g
CC_FLAGS := -g -O3
CFLAGS   := -O2 -g -Wall
CFLAGS_1 := -ansi -std=gnu99
CFLAGS_2 := -ansi -pedantic -std=gnu99 -Wall -W
CFLAGS_3 := -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -g -O2

# Valgrind
VAL      := valgrind --track-origins=yes -v
VALFULL  := valgrind --track-origins=yes --leak-check=full -v
VALMEM   := valgrind --track-origins=yes --tool=memcheck --leak-check=full --show-leak-kinds=all --show-reachable=yes --num-callers=20 --track-fds=yes -v
VALMASS  := valgrind --tool=massif prog

# c files
# SOURCES  := $(CFILES) # $(CFILES2) all CFILES
CPPFILES   := $(wildcard src/*.cpp)


# o files
OBJECTS  := $(SOURCES:.cpp=.o)
OBJECTS_A:= $(SOURCES:.cpp=_a.o)

# lib
LIB      := -pthread
# LIB    := -pthread -larmadillo
# LIB    := -pthread -lmongoclient -L lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt

# include
INC      := -Iinclude
INC_EIGEN:= -I/usr/include/eigen3


# DIRS
TOPDIR   :=	$(shell pwd)
SRCDIR   := src
BUILDDIR := build
BIN_DIR  := bin
BIN      := run_segment
SRCEXT   := cpp
ODIR     := build
TESTD    := test

# Executable:
EXEC     := run_segment
EXEF     := $(wildcard /usr/local/bin/run_segment*)


# ------------------------------------------------------------------------------
# For molecule reading: ONEMOL, DEFAULT, DCDREAD, ALLATOM
# Macros
MACRO = -D
INFO  = -DINFO # optional
DCD = -DDCDREAD
# Original: unmodified, like a subset of frames
USEDCD = -DDCDREAD -DDCD_WRITE_B -DDCD_WRITE -DDCD_WRITE_UNMOD -DDCD_WRITE_E
USEDCDM= -DDCDREAD -DDCD_WRITE_B -DDCD_WRITE -DDCD_WRITE_E
# USEDCDR= -DDCDREAD -DDCD_WRITE_B -DDCD_WRITE -DDCD_WRITE_ROT -DDCD_WRITE_E

# SM: tension/costheta,contacts-MAPSM
TEN   = $(DCD) -DTENSION_COSTHETA
# CHI   = $(DCD) -DCHIVALUE
CHI   = $(DCD) -DCHI_BEGIN -DCHI_MID -DCHI_END

# Other. indices,
INDICES = -DONEMOL -DINDICES
INDICES_AA = -DALLATOM -DINDICES_AA
TOP_CHARMM = -DALLATOM -DTOP_CHARMM

# 2 mol
DEFAULT= -DDEFAULT

# 1 mol
ONEMOL = -DONEMOL
TOPO   = -DONEMOL -DTOPOLOGY
TOPOPF = -DONEMOL -DTOPOPF
MTTOP  = -DONEMOL -DMTTOP

# contacts: MT uses LonN,LatE etc. SM uses chain_ref[0].contacts (intra)
# MT: with latlon
# PF: with lononly_proto
# single molecule (no latlon)
# CONFILE    = -DCONTACTFILE
CONFILE    = -DCONFILE
CONTACTBYRESFILE1 = -DCONTACTSBYRESFILE1
CONMAP_MT  = -DCONTACTMAP_MT1 -DCONTACTMAP_MT2 -DCONTACTMAP_MT3
CONMAP_SM  = -DCONTACTMAP_SM1 -DCONTACTMAP_SM2 -DCONTACTMAP_SM3
ANG3       = -DANGLE3CENTROID_B -DANGLE3CENTROID_M -DANGLE3CENTROID_E
CENMT      = -DCENTROIDMOVEMENT_B -DCENTROIDMOVEMENT_M -DCENTROIDMOVEMENT_E
CURVATURE3 = -DCURVATURE_B -DCURVATURE_M -DCURV3 -DCURVATURE_E
# CONANGCURV = -DCONANGCURV_B -DCONANGCURV_M -DCONANGCURV_E
DIMERMAPDEF= -DDIMERMAP -DDIMERMAP_1A -DDIMERMAP_1B -DDIMERMAP_2 -DDIMERMAP_3
# DIMERCURV  = -DDMCURV_B -DDMCURV_M -DDMCURV_E # requires DIMERMAPDEF
CURVATURE4 = -DCURVATURE_B -DCURVATURE_M -DCURV4 -DCURVATURE_E
MTPF_SECS  = -DMTPF_B -DMTPF_M -DMTPF_E
# MTDIMER_S  = -DMT_DIMER_CONTACT_B -DMT_DIMER_CONTACT_M -DMT_DIMER_CONTACT_E

# MT: contacts, centroid move, angle, con|angle|cenmov
# MT, PF & SM
MAPPF       = $(DCD) -DLONONLY_PROTO $(CONMAP_MT) $(CONFILE)
MAPSM       = $(DCD) $(CONMAP_SM) $(CONTACTBYRESFILE1) # $(CONFILE)
MAPMT       = $(DCD) -DLATLON $(CONMAP_MT) $(CONFILE)
CENMOV      = $(DCD) -DLATLON $(CENMT)
ANGLEC_MT   = $(DCD) -DLATLON $(ANG3)
CON_ANG_CURV= $(DCD) -DLATLON $(ANG3) $(CONMAP_MT) $(CURVATURE3) $(CONFILE)
DIMERMAP    = $(DCD) -DLATLON $(CONMAP_MT) $(DIMERMAPDEF) # note. no CONFILE
MTPF        = $(DCD) -DLATLON $(MTPF_SECS)
# MTPF = $(DCD) -DLONONLY_PROTO $(MTPF_SECTIONS) # ? maybe
# MTDIMER     = $(DCD) -DLATLON $(MTDIMER_S)


# PF: cm,curva,angles,contacts-MAPPF),proto curvature & centroid movement
ANGLEC_PF   = $(DCD) -DLONONLY_PROTO $(ANG3)
# proto angles (by resids, 14 257) ---> suspect!!! vvv
PFANGLE   = $(DCD) -DLONONLY_PROTO -DPFANGLE


# ----- deprecated ---------------------------------------------------------
# DCDPROTO= -DDCDREAD -DPROTO
# INERTIA = -DONEMOL -DINERTIA

# PROTO = -DONEMOL -DPROTO
# DCDTOP= -DDCDREAD -DTOPOLOGY
# DCDINDEX= -DDCDREAD -DINDICES

# Deprecated
# -------
# suspect: CONTACT, CONTACT_PERSIST
# DCDCON = -DDCDREAD -DCONTACT -DCONTACTPERSIST -DTENSION_COSTHETA
# CONTACT= -DDEFAULT -DCONTACT -DCONTACTPERSIST -DTENSION_COSTHETA
# =======
# MTCON2, MTCON_OC,MTCON1
# =======
# DCDMT = -DDCDREAD -DMTCON1 -DMTCON2
# ------> fails - was using "compare_contacts",now uses "evaluate_original_contacts_now
# MTCONTACT = -DDEFAULT -DMTCON1 -DMTCON2
# DCD_MTCON_OC = -DDCDREAD -DMTCON1 -DMTCON_OC


# ------------------------------------------------------------------------------
# Main target
$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -lm -o $(EXEC)

# To remove generated files
clean:
	rm -f $(OBJECTS) $(OBJECTS_A)
# rm -f $(EXEC) $(OBJECTS)

go:
	./$(EXEC)

foo:
	@echo echoing foo
	@echo $(CPPFILES)
	@echo $(SOURCES)
	@echo $(OBJECTS)
	@echo $(OBJECTS_A)

# ---------------------------------------------------------------------
# Examples & Testing
# ---------------------------------------------------------------------
def: $(OBJECTS)
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DEFAULT) -o test/$(EXEC)_def
	cd test && $(EXEC)_def proto1.pdb proto2501.pdb 14 2
main:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DEFAULT) -o test/$(EXEC)_def
	cd test && ./$(EXEC)_def proto1.pdb proto2501.pdb 14 2

# costheta, tension: verified November 23, 2015.
costension_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TEN) -DNDEBUG -o bin/$(EXEC)_dcd_costension_1
costension:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TEN) -o test/$(EXEC)_dcd_costension_1
	cd test && ./$(EXEC)_dcd_costension_1 nbd_ref.pdb nbd_2KHO_D30_pull.dcd 10 1 0
# cd test && ./$(EXEC)_dcd_contact nbd_2KHO.ref.pdb nbd_2KHO_D30_pull.dcd 5 1 0

indices_aa_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(INDICES_AA) -DNDEBUG -o bin/$(EXEC)_index_aa
indices_aa:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(INDICES_AA) -o test/$(EXEC)_index_aa
	cd test && ./$(EXEC)_index_aa 4EZW_no-anisou.pdb 19 0
indices_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(INDICES) -DNDEBUG -o bin/$(EXEC)_index
indices:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(INDICES) -o test/$(EXEC)_index
	cd test && ./$(EXEC)_index mtbig.pdb 209 1 37
topcharmm_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOP_CHARMM) -DNDEBUG -o bin/$(EXEC)_topcharmm
topcharmm:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOP_CHARMM) -o test/$(EXEC)_topcharmm
	cd test && ./$(EXEC)_topcharmm 4EZW_no-anisou.pdb 19 0
# cd test && ./$(EXEC)_topcharmm 2KHO.pdb 1 0
# cd test && ./$(EXEC)_topcharmm 2khosbd_caonly.pdb 1 0

# PF curva/cenmov, use DCD, USEDCDM
# ---run_segment_dcd_curva_pfm
# ---analyses
# curvature_global.dat
# curvature_local.dat
# distance_proto_centroids.dat
# radius_curvature_global.dat
# radius_curvature_local.dat
curva_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(CURVATURE4) -DNDEBUG -o bin/$(EXEC)_dcd_curva4_pf
curva:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) $(CURVATURE4) -o test/$(EXEC)_dcd_curva4_pf
	cd test && ./$(EXEC)_dcd_curva4_pf proto_interim_f1_r-1-5641.pdb proto_TUB_multiEh_k23_N1_pull.dcd 100 14 2 14 257
curva_mod_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCDM) $(CURVATURE4) -DNDEBUG -o bin/$(EXEC)_dcd_curva4_pfm
curva_mod:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCDM) $(CURVATURE4) -o test/$(EXEC)_dcd_curva4_pfm
	cd test && time ./$(EXEC)_dcd_curva4_pfm kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 50 13900 900
# cd test && time ./$(EXEC)_dcd_curva_pfm kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 3000 14000 1543
curva_val:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCDM) $(CURVATURE4) -o test/$(EXEC)_dcd_curva4_pfm
	cd test && $(VALFULL) ./$(EXEC)_dcd_curva4_pfm kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 50 15 3 14 257
curva_mod3_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCDM) $(CURVATURE3) -DNDEBUG -o bin/$(EXEC)_dcd_curva3_pfm
curva_mod3:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCDM) $(CURVATURE3) -o test/$(EXEC)_dcd_curva3_pfm
	cd test && time ./$(EXEC)_dcd_curva3_pfm kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 50 13900 900

# PF angle: verified, November 23, 2015.
pfangle_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(PFANGLE) -DNDEBUG -o bin/$(EXEC)_dcd_angle_pf
pfangle:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(PFANGLE) -o test/$(EXEC)_dcd_angle_pf
	cd test && ./$(EXEC)_dcd_angle_pf kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 150 15 3 14 257

# MT angle: verified, November 23, 2015.
# ---run_segment_dcd_anglec_mt/pf
# ---testing----
mtpf_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MTPF) -DNDEBUG -o bin/$(EXEC)_dcd_mtpf
mtpf:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MTPF) -o test/$(EXEC)_dcd_mtpf
	cd test && time ./$(EXEC)_dcd_mtpf mt.ref.pdb mt_d100_indent.dcd 209 1 0 500 25
# mtpf1: # not ready
# 	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MTPF) -o test/$(EXEC)_dcd_mtpf_pf
# 	cd test && time ./$(EXEC)_dcd_mtpf_pf kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 0 20000 1500
# 	# cd test && time ./$(EXEC)_dcd_mtpf mt.ref.pdb mt_d100_indent.dcd 209 1 0 500 25
# ---analyses---
# mt_angles_ew.dat,mt_angles_ns.dat
anglec_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(ANGLEC_MT) -DNDEBUG -o bin/$(EXEC)_dcd_anglec_mt
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(ANGLEC_PF) -DNDEBUG -o bin/$(EXEC)_dcd_anglec_pf
anglec_mt:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(ANGLEC_MT) -o test/$(EXEC)_dcd_anglec_mt
	cd test && time ./$(EXEC)_dcd_anglec_mt mt.ref.pdb mt_d100_indent.dcd 209 1 0 500 25
anglec_pf: # to protofilaments
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(ANGLEC_PF) -o test/$(EXEC)_dcd_anglec_pf
	cd test && time ./$(EXEC)_dcd_anglec_pf kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 100 800 60
# kinesin13_123.ref.pdb
# cd test && time ./$(EXEC)_dcd_angle_mt mt.ref.pdb mt_d100_indent.dcd 209 1 0 500 25
# cd test && ./$(EXEC)_dcd_angle_mt mtonly_seamup.ref.pdb mtonly_seamup_d1_indent.dcd 85 209 1

# MT only: verified, November 23, 2015.
cenmov_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(CENMOV) -DNDEBUG -o bin/$(EXEC)_dcd_cenmov_mt
cenmov:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(CENMOV) -o test/$(EXEC)_dcd_cenmov_mt
	cd test && time ./$(EXEC)_dcd_cenmov_mt mt.ref.pdb mt_d100_indent.dcd 209 1 200 1000 80
# cd test && ./$(EXEC)_dcd_cenmov_mt mtonly_seamup.ref.pdb mtonly_seamup_d1_indent.dcd 85 209 1


# contacts: verified all 3, November 23, 2015.
# ---run_segment_dcd_contactmap_| mt | pf | 1 |.
# ---analyses---CONTACTFILE
# contacts_e.dat
# contacts_n.dat
# contacts_s.dat
# contacts_w.dat
contactmapM_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MAPMT) -DNDEBUG -o bin/$(EXEC)_dcd_contactmap_mt
contactmapP_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MAPPF) -DNDEBUG -o bin/$(EXEC)_dcd_contactmap_pf
contactmapN_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MAPSM) -DNDEBUG -o bin/$(EXEC)_dcd_contactmap_1
contactmapM:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MAPMT) -o test/$(EXEC)_dcd_contactmap_mt
	cd test && time ./$(EXEC)_dcd_contactmap_mt kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 100 800 60
contactmapdimer:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MAPMT) -o test/$(EXEC)_dcd_contactmap_mt
	cd test && time ./$(EXEC)_dcd_contactmap_mt mt.ref.pdb mt_d100_indent.dcd 209 1 0 500 125
contactmapP:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MAPPF) -o test/$(EXEC)_dcd_contactmap_pf
	cd test && time ./$(EXEC)_dcd_contactmap_pf kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 100 13000 900
contactmapN:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MAPSM) -o test/$(EXEC)_dcd_contactmap_1
	cd test && ./$(EXEC)_dcd_contactmap_1 2khosbd_16-4.ref.pdb 2khosbd_16-4.ref.dcd 1 0 100 800 40
# cd test && time ./$(EXEC)_dcd_contactmap_1 kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 100 800 60

# # combined: angles & contacts.
# dimermap_release:
# 	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(CON_ANG_CURV) -DNDEBUG -o bin/$(EXEC)_dcd_dimermap_mt
# dimermap:
# 	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(CON_ANG_CURV) -o test/$(EXEC)_dcd_dimermap_mt
# 	cd test && time ./$(EXEC)_dcd_dimermap_mt mt.ref.pdb mt_d100_indent.dcd 209 1 0 500 125
# combined: angles & contacts.
dimermap_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DIMERMAP) -DNDEBUG -o bin/$(EXEC)_dcd_dimermap_mt
dimermap:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DIMERMAP) -o test/$(EXEC)_dcd_dimermap_mt
	cd test && time ./$(EXEC)_dcd_dimermap_mt mt.ref.pdb mt_d100_indent.dcd 209 1 0 800 125
# cd test && $(VALFULL) ./$(EXEC)_dcd_dimermap_mt mt.ref.pdb mt_d100_indent.dcd 209 1 0 800 125

dimermap_reverse:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DIMERMAP) -o test/$(EXEC)_dcd_dimermap_mt
	cd test && time ./$(EXEC)_dcd_dimermap_mt mt.ref545.pdb mt_retract545.dcd 209 1 0 500 30
# cd test && time ./$(EXEC)_dcd_dimermap_mt mt_lat3.pdb mt_retract.dcd 209 1 0 500 125
dimermap_reversev:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DIMERMAP) -o test/$(EXEC)_dcd_dimermap_mt
	cd test && time $(VALFULL) ./$(EXEC)_dcd_dimermap_mt mt_lat3.pdb mt_retract.dcd 209 1 0 500 125
dimermap302:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DIMERMAP) -o test/$(EXEC)_dcd_dimermap_mt
	cd test && time ./$(EXEC)_dcd_dimermap_mt mt_lat3_02.pdb mt_lat3_02.dcd 209 1 0 500 125


# chi
chi_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(CHI) -DNDEBUG -o bin/$(EXEC)_dcd_chi_1
chi:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(CHI) -o test/$(EXEC)_dcd_chi_1
	cd test && ./$(EXEC)_dcd_chi_1 2khosbd_16-4.ref.pdb 2khosbd_16-4.ref.dcd 1 0 5 10005 5
# cd test && ./$(EXEC)_dcd_chi_1 2khosbd_16-4.ref.pdb 2khosbd_16-4.ref.dcd 1 0 5 10005 5 0 220

# use_dcd:
write_dcd_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCD) -DNDEBUG -o bin/$(EXEC)_dcd_write
write_dcd:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCD) -o test/$(EXEC)_dcd_write
	cd test && time ./$(EXEC)_dcd_write kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 20 13000 50
# mod_dcd_release:
# 	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCDM) -DNDEBUG -o bin/$(EXEC)_dcd_rotxy
# mod_dcd:
# 	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(USEDCDM) -o test/$(EXEC)_dcd_write_rotxy
# 	cd test && ./$(EXEC)_dcd_write_rotxy kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 150 15 3
# # cd test && ./$(EXEC)_dcd_writemod nbd_ref.pdb nbd_2KHO_D30_pull.dcd 150 1 0

# topology
topology_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPO) -DNDEBUG -o bin/$(EXEC)_top
topology:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPO) -DSBDPEP -o test/$(EXEC)_top_testpep
	cd test && ./$(EXEC)_top_testpep 1_DH_ca.pdb 2 0
topology_release_pf:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOPF) -DNDEBUG -o bin/$(EXEC)_top-pf
topology_pf:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPOPF) -o test/$(EXEC)_top-pf
	cd test && ./$(EXEC)_top-pf kinesin13_123.ref.pdb 15 3
topology-nopep:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPO) -DSBDNOPEP -o test/$(EXEC)_top_nopep
	cd test && ./$(EXEC)_top_nopep pdb2khosbd_CA.ent 1 0
topology-sbdchainA:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(TOPO) -DSBDNOPEPchainA -o test/$(EXEC)_top_sbdchainA
topology_mt_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MTTOP) -DNDEBUG -o bin/$(EXEC)_top_mt
topology_mt:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(MTTOP) -o test/$(EXEC)_top_mt
	cd test && ./$(EXEC)_top_mt mtonly_seamup.ref.pdb 209 1
dcd_release:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -DNDEBUG -o bin/$(EXEC)_dcd
dcd:
	$(CXX) $(CPPFILES) $(CF) $(INC) $(LIB) $(DCD) -o test/$(EXEC)_dcd
	cd test && time ./$(EXEC)_dcd kinesin13_123.ref.pdb kinesin13_123_D1_pull.dcd 15 3 3000 12503 10

# -----------------------------------------------------------------------------
# Make all.
all: \
	dimermap_release \
	indices_release \
	topology_release \
	write_dcd_release \
	costension_release \
	curva_mod_release \
	curva_mod3_release \
	anglec_release \
	cenmov_release \
	contactmapM_release \
	contactmapP_release \
	contactmapN_release \
	dcd_release \
	chi_release


# curva_release \ pf(bad) not pfm(good)
# pfangle_release \ # uses the 14,257 indices; suspect.
# anglemt_release \

# install:
# 	echo $(EXEF)
# 	$(shell rm) $(EXEF)

# verified recently?
# contact_release: needs fail safe.
# proto_ : need the atom indices

# suspect.
# inertia_release
# proto_release \

# Deprecated.
# contact_release \
# contact_dcd_release \

# -----------------------------------------------------------------------------
# To obtain object files (needs to be at the end??)
# %.o: %.cpp
# 	$(CXX) $(CC_FLAGS) $(INC) $(LIB) -c $< -o $@
