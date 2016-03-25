
CC = gcc
STRIP = strip
RM = rm
CFLAGS = -O3 -std=gnu99 -D__RUN_IN_LINUX__ -march=native
OBJECTS = dynamicMem.o random.o GA.o timer.o LABS.o TabuSearch.o
UNAME := $(shell uname -s)

ifeq ($(SSE4.2),true)
        CFLAGS:=$(CFLAGS) -msse4.2
endif

ifeq ($(MAIN),lssRRts)
        CFLAGS:=$(CFLAGS) -DlssRRts
endif

ifeq ($(MAIN),lssRRts0)
        CFLAGS:=$(CFLAGS) -DlssRRts0
endif

EXE=$(MAIN)
ifeq ($(M32),true)
        CFLAGS:=$(CFLAGS) -m32
	EXE:=$(MAIN).32bit
endif

ifeq ($(M64),true)
        CFLAGS:=$(CFLAGS) -m64
endif

LINK = gcc -o $(EXE) $(OBJECTS) -g0 -lpthread libs/dcmt0.6.1/lib/libdcmt.a $(CFLAGS)

all: LIB
	$(RM) -f *.o lssMAts *~
	make lssMAts$(ARCH) MAIN=lssMAts CFLAGS="$(CFLAGS) -g0 -Ofast -DNDEBUG -DlssMAts" 
	$(STRIP) lssMAts
	$(RM) -f *.o lssRRts *~
	make lssRRts MAIN=lssRRts CFLAGS="$(CFLAGS) -g0 -Ofast -DNDEBUG -DlssRRts" 
	$(STRIP) lssRRts
	#$(RM) -f *.o lssRRts0 *~
	#make lssRRts0 MAIN=lssRRts0 CFLAGS="$(CFLAGS) -g0 -Ofast -DNDEBUG -DlssRRts0" 
	#$(STRIP) lssRRts0

.c.o: *.c *.h
	$(CC) -c $(CFLAGS) $<

$(MAIN): $(OBJECTS)
	$(LINK)

$(MAIN).32bit: $(OBJECTS)
	$(LINK)
	$(STRIP) $(EXE)

debug: LIB
	$(RM) -f *.o lssMAts *~
	make lssMAts$(ARCH) MAIN=lssMAts CFLAGS="$(CFLAGS) -g3 -DlssMAts"
	$(RM) -f *.o lssRRts *~
	make lssRRts$(ARCH) MAIN=lssRRts CFLAGS="$(CFLAGS) -g3 -DlssRRts"
	$(RM) -f *.o lssRRts0 *~
	make lssRRts0$(ARCH) MAIN=lssRRts0 CFLAGS="$(CFLAGS) -g3 -DlssRRts0"

LIB:
	make M32=$(M32) M64=$(M64) -C libs/dcmt0.6.1/lib/

clean:
	make clean -C libs/dcmt0.6.1/lib/
	$(RM) -f *.o lssMAts lssRRts lssRRts0 *~
