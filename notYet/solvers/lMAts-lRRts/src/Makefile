
CC=gcc 
CFLAGS=-g0 -std=gnu99 -O3 -DNDEBUG -D__RUN_IN_LINUX__ -march=native
LDFLAGS=-g0 -lpthread -march=native
OBJECTS=dynamicMem.o random.o GA.o timer.o LABS.o TabuSearch.o libs/dcmt0.6.1/lib/libdcmt.a
UNAME := $(shell uname -s)
MAIN=lMAts
STRIP = strip

all: libs/dcmt0.6.1/lib/libdcmt.a
	$(RM) -f *.o lMAts *~
	make lMAts MAIN=lMAts CFLAGS="$(CFLAGS) -g0 -Ofast -DNDEBUG -DlMAts"
	$(RM) -f *.o lRRts *~
	make lRRts MAIN=lRRts CFLAGS="$(CFLAGS) -g0 -Ofast -DNDEBUG -DlRRts"

.c.o:
	$(CC) -c $(CFLAGS) $<

libs/dcmt0.6.1/lib/libdcmt.a:
	make -C libs/dcmt0.6.1/lib/

$(MAIN): $(OBJECTS)
	gcc -o $(MAIN) $(OBJECTS) $(LDFLAGS)
	strip $(MAIN)

clean:
	rm -f *.o *~ lMAts lRRts
	rm -f $(MAIN)
	make clean -C libs/dcmt0.6.1/lib/

