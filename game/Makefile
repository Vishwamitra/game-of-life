EXECUTABLES=game
EXPENSIVE_JUNK += $(EXECUTABLES)
SRC = game.c

JUNK +=

CFLAGS += -O3 -Wall -W --std=c11 -lm
CXXFLAGS += -O3 -Wall -W --std=c++11 -lm -Wno-cast-function-type
OMP_CFLAGS = $(CFLAGS) -fopenmp
MPI_CFLAGS = $(CXXFLAGS) -lmpi

all: $(EXECUTABLES)

clean:
		rm -rf $(JUNK)