******************************************************
The Low Autocorrelation Binary Sequence (LABS) Problem
******************************************************

Compile:
=======

make

It works on GNU/Linux, MAC, and other systems.

GCC 4.4 or higher is required.

For new processors use the following parameter: SSE4.2=true

For large sequence sizes (>320) use the following parameter: MAX_L=value
The value must be a multiple of 64.

Examples:

 make
 make SSE4.2=true
 make MAX_L=1024

For debug version use the following target: debug

 make debug

Debug version has additionaly flagis: -trace and -walk

To run the program for LABS problem instance of L=43 (on Linux) write:

 ./lssOrel 43 -v


Some examples on how to use and some outputs are in testRuns directory


