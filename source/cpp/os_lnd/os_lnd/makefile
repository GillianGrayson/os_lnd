CC = icc
CFLAGS = -O2 -mkl=sequential -qopenmp -std=c++14

PETSC_DIR=/home/ivanchen/yusipov/os_lnd/source/cpp/os_lnd/os_lnd
PETSC_ARCH=.
SLEPC_DIR=/home/ivanchen/yusipov/os_lnd/source/cpp/os_lnd/os_lnd

RM  = rm -f
CP  = cp
MV  = mv

os_lnd:
	$(CC) $(CFLAGS) -I./third_party/eigen-git-mirror -I ./boost_1_71_0 -I ./third_party/fmt/include -I ./third_party/gcem/include -I ./third_party/inih -I ./third_party/odeint-v2/include -I ./third_party/spdlog/include -I ./third_party/filesystem/include libslepc.so lib/libpetsc.so -I$(PETSC_DIR)/include -I$(SLEPC_DIR)/include main.cpp -o os_lnd.o
	
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common