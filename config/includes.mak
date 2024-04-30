CXX:=g++-13
CC:=gcc-13
OS:=OSX
CFLAGS:= -march=native -O3 -fopenmp -DGSL -DKNUTH
CFLAGS_EXHIB:= -march=native -O3 -fopenmp -DGSL -DHAVE_INLINE

IBASIC:=-I/usr/local/include/
LINK:=-L/usr/local/lib/
LIBS_EXHIB:=-lm -lgsl -lgslcblas

# directories

LD_EXHIB_SRC:= Lie_Detector_exhibition/

LD_EXHIB_DIR:= $(addsuffix objs/, $(LD_EXHIB_SRC))

TEST_SRC:=tests/

# limited include directory, for demonstration
IDIR_EXHIB:=$(IBASIC) -I$(LD_EXHIB_SRC)
