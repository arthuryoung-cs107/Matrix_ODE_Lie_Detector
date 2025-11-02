CXX:=g++-15
CC:=gcc-15
OS:=OSX
CFLAGS_EXHIB:= -march=native -O3 -fopenmp -DGSL -DHAVE_INLINE

IBASIC:=-isystem/usr/local/include/
LINK:=-L/usr/local/lib/
LIBS_EXHIB:=-lm -lgsl -lgslcblas
# LIBS_EXHIB:=-lm -lgsl -lcblas

# directories

LD_EXHIB_SRC:= Lie_Detector_exhibition/

LD_EXHIB_DIR:= $(addsuffix objs/, $(LD_EXHIB_SRC))

TEST_SRC:=tests/

# limited include directory, for demonstration
IDIR_EXHIB:=$(IBASIC) -I$(LD_EXHIB_SRC)
