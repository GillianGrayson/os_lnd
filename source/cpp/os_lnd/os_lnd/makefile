CC = icc
CFLAGS = -O2 -mkl=sequential -qopenmp -std=c++14

PETSC_DIR=/home/ivanchen/yusipov/eigs/petsc-3.12.4
PETSC_ARCH=arch-linux2-c-opt
SLEPC_DIR=/home/ivanchen/yusipov/eigs/slepc-3.12.2

RM  = rm -f
CP  = cp
MV  = mv

os_lnd:
	$(CC) $(CFLAGS) -I./third_party/eigen -I ./boost_1_71_0 -I ./third_party/fmt/include -I ./third_party/gcem/include -I ./third_party/inih -I ./third_party/odeint-v2/include -I ./third_party/spdlog/include -I ./third_party/filesystem/include -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include  -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include main.cpp -o os_lnd.o
	
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common