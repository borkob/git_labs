CFLAGS = -Wall -O3 -I../include -L../lib
CC = gcc

all:example1 example2 example3 new_example1 new_example2 new_example3

example1:   example1.c ../include/dc.h ../lib/libdcmt.a
	$(CC) $(CFLAGS) -o example1 example1.c -ldcmt

example2:   example2.c ../include/dc.h ../lib/libdcmt.a
	$(CC) $(CFLAGS) -o example2 example2.c -ldcmt

example3:   example3.c ../include/dc.h ../lib/libdcmt.a
	$(CC) $(CFLAGS) -o example3 example3.c -ldcmt

new_example1: new_example1.c ../include/dc.h ../lib/libdcmt.a
	$(CC) $(CFLAGS) -o new_example1 new_example1.c -ldcmt

new_example2: new_example2.c ../include/dc.h ../lib/libdcmt.a
	$(CC) $(CFLAGS) -o new_example2 new_example2.c -ldcmt

new_example3:  new_example3.c ../include/dc.h ../lib/libdcmt.a
	$(CC) $(CFLAGS) -o new_example3 new_example3.c -ldcmt

clean :
	/bin/rm -f *.o example1 example2 example3 new_example1 new_example2 \
		       new_example3


