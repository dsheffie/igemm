OBJ = gemm.o igemm.o
CXX = riscv32-unknown-elf-g++
CC = riscv32-unknown-elf-gcc
EXE = gemm
OPT = -O3
CXXFLAGS = -std=c++11 -g $(OPT)
CFLAGS = -std=c11 -g $(OPT)
DEP = $(OBJ:.o=.d)

.PHONY: all clean

all: $(EXE)

$(EXE) : $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LIBS) -o $(EXE) -specs=htif.specs

%.o: %.cc
	$(CXX) -MMD $(CXXFLAGS) -c $<

%.o: %.c
	$(CC) -MMD $(CFLAGS) -c $< 

-include $(DEP)

clean:
	rm -f $(EXE) $(OBJ) $(DEP) *.csv *~
