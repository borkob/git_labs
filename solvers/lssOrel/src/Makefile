CXX=g++
UNAME := $(shell uname -s)
DATE:=$(shell date '+%d_%m_%G')
RAND:=$(shell /bin/bash -c "echo $$RANDOM")
DIR:=$(DATE)_$(RAND)
VERSION=004
NAME=lssOrel
BIN=$(NAME)

ifndef MAX_L
MAX_L=320
endif

CXXFLAGS=-Wall -W -std=c++0x -DMAX_L=$(MAX_L) -DNAME='\"$(NAME)\"' -DVERSION='\"$(VERSION)\"' -march=native

ifeq ($(SSE4.2),true)
	CXXFLAGS:=$(CXXFLAGS) -msse4.2
endif

ifeq ($(M32),true)
	CXXFLAGS:=$(CXXFLAGS) -m32
	BIN:=$(NAME).32bit
endif

ifeq ($(M64),true)
	CXXFLAGS:=$(CXXFLAGS) -m64
endif

OBJ=sequence.o saw.o main.o

all:
	rm -f *.o $(BIN)
	make $(BIN) CXXFLAGS="$(CXXFLAGS) -DNDEBUG -g0 -O3"
	strip $(BIN)

debug:
	make $(BIN) CXXFLAGS="$(CXXFLAGS) -g3"

profile:
	make $(BIN) CXXFLAGS="$(CXXFLAGS) -DNDEBUG -pg -O3"

profile2:
	make $(BIN) CXXFLAGS="$(CXXFLAGS) -DNDEBUG -fprofile-arcs -ftest-coverage"

$(BIN): $(OBJ) *.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(BIN)

$(BIN)B: $(OBJ) *.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(BIN)B

%.o : %.cpp *.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o $(BIN) lssOrel *.gcda *.gcno *.gcov
