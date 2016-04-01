############################################## COMPILE ######################################################

$ make clean
$ make

# for random strategy:
$ make STRATEGY=random

# for debuging (trace option is added)
$ make debug

############################################################################################################

$ bin/lssMAts

usage: bin/lssMAts <binary seq length> <random seed> <max time (secs)> <no. of threads> <valueTarget>

note1: The last argument is optional. If no <valueTarget> is specified,
        the "best known value", stored internally, will be accessed.
note2: This program stops either
        if runtime   >= <max time (secs)> or
        if valueBest <= <valueTarget>

Copyright 2012
*  José E. Gallardo, Carlos Cotta, and Antonio J. Fernández
*  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
*  Applied Soft Computing. 9(4): 1252-1262 (2009).
*  ---------------------------------------------------------------------
*  Original source code was modified /instrumented for testing by Borko Boskovic.

###############################################################################################################

$ bin/lssRRts 
usage: bin/lssRRts <binary seq length> <random seed> <max time (secs)> <no. of threads> <valueTarget>

note1: The last argument is optional. If no <valueTarget> is specified,
        the "best known value", stored internally, will be accessed.
note2: This program stops either
        if runtime   >= <max time (secs)> or
        if valueBest <= <valueTarget>

Copyright 2012
*  José E. Gallardo, Carlos Cotta, and Antonio J. Fernández
*  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
*  Applied Soft Computing. 9(4): 1252-1262 (2009).
*  ---------------------------------------------------------------------
*  Original source code was modified /instrumented for testing by Borko Boskovic.
