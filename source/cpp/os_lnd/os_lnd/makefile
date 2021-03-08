CC = icc
CFLAGS = -O2 -mkl=sequential -qopenmp -std=c++14

RM  = rm -f
CP  = cp
MV  = mv

os_lnd:
	$(CC) $(CFLAGS) -I./third_party/eigen -I ./boost_1_71_0 -I ./third_party/fmt/include -I ./third_party/gcem/include -I ./third_party/inih -I ./third_party/odeint-v2/include -I ./third_party/spdlog/include -I ./third_party/filesystem/include main.cpp -o os_lnd.o