/*
  lssOrel, a solver for LABS problem
  Copyright (C) 2014 Borko Boršković, Franc Brglez and Janez Brest

  lssOrel is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  lssOrel is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "saw.h"

#include <iostream>
#include <sstream>
#include <libgen.h>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

const unsigned int MAX_BEST_L = 305;
const unsigned bestT[MAX_BEST_L] = {
       0,    0,    1,    1,    2,    2,    7,    3,    8,   12,   13,    5,
      10,    6,   19,   15,   24,   32,   25,   29,   26,   26,   39,   47,
      36,   36,   45,   37,   50,   62,   59,   67,   64,   64,   65,   73,
      82,   86,   87,   99,  108,  108,  101,  109,  122,  118,  131,  135,
     140,  136,  153,  153,  166,  170,  175,  171,  192,  188,  197,  205,
     218,  226,  235,  207,  208,  240,  257,  241,  250,  274,  295,  275,
     300,  308,  341,  329,  334,  358,  347,  339,  352,  372,  377,  377,
     430,  414,  439,  431,  448,  432,  453,  477,  498,  486,  499,  479,
     520,  536,  545,  577,  578,  578,  567,  555,  612,  620,  701,  677,
     702,  662,  723,  687,  788,  752,  817,  745,  814,  786,  847,  835,
     872,  844,  885,  893,  922,  846,  875,  887,  932,  920,  945,  913,
    1014, 1010, 1063, 1027, 1076, 1052, 1117, 1133, 1178, 1126, 1235, 1191,
    1248, 1208, 1273, 1265, 1298, 1218, 1279, 1275, 1388, 1340, 1429, 1437,
    1438, 1366, 1467, 1439, 1504, 1512, 1585, 1529, 1594, 1474, 1591, 1563,
    1620, 1532, 1693, 1677, 1674, 1606, 1699, 1719, 1780, 1808, 1929, 1897,
    1898, 1898, 2015, 1995, 2040, 2028, 2069, 1973, 1970, 1966, 2123, 2191,
    2304, 2272, 2337, 2281, 2374, 2218, 2343, 2275, 2412, 2460, 2541, 2421,
    2542, 2662, 2723, 2695, 2720, 2664, 2761, 2801, 2878, 2698, 2799, 2831,
    2968, 3036, 3173, 3189, 3322, 3206, 3211, 3215, 3392, 3416, 3569, 3409,
    3566, 3474, 3687, 3587, 3752, 3692, 3821, 3757, 3674, 3590, 3651, 3711,
    3948, 3992, 4073, 4073, 4150, 4098, 4223, 4291,    0, 4280, 4341, 4165,
    4386, 4382, 4587, 4463, 4584, 4472, 4705, 4705,    0, 4790, 4887, 4803,
    4892, 4948,    0, 5037, 5130, 4950, 5203, 5243,    0,    0,    0,    0,
       0,    0,    0,    0,    0, 5564,    0, 5697,    0, 5790,    0,    0,
       0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
       0,    0,    0, 6455, 6568
};

void help(const char argv[]){
    const unsigned int w = 22;
    cout<<left<<"Usage: "<<argv<<" Lodd [options]"<<endl<<endl;
    cout<<setw(w)<<"Lodd";
    cout<<"... odd number, denoting the length of the binary string"<<endl;
    cout<<endl<<"Options:"<<endl;
    cout<<setw(w)<<"-runtimeLmt seconds";
    cout<<"... stop after this number of seconds; default = 300."<<endl;
    cout<<setw(w+4)<<""<<"(this can be overruled if cntProbeLmt below is enabled)"<<endl;
    cout<<setw(w)<<"-valueTarget energy";
    cout<<"... stop when this level of energy is reached;"<<endl;
    cout<<setw(w+4)<<""<<"default = best known energy value."<<endl;
    cout<<setw(w)<<"-laevus number";
    cout<<"... functionFactor as well as functionParameter."<<endl;
    cout<<setw(w+4)<<""<<"The laevus value represents the number of 0 bits at the left end"<<endl;
    cout<<setw(w+4)<<""<<"of full binary string; default = 0."<<endl;
    cout<<setw(w+4)<<""<<"NOTE: since this solver runs under the assumption of skew-symmetry,"<<endl;
    cout<<setw(w+10)<<""<<"the effective dimension of the search space is reduced to"<<endl;
    cout<<setw(w+16)<<""<<"nDim = (Lodd + 1)/2 - laevus"<<endl;
    cout<<setw(w+10)<<""<<"which now controls the values of solver parameters."<<endl;
    cout<<setw(w)<<"-walkSegmCoef val"<<"... solverFactor, a coefficient that controls walk segment length;"<<endl;
    cout<<setw(w+4)<<""<<"default = 8 implies walk segment length value of 8*(Lodd + 1)/2 - laevus"<<endl;
    cout<<setw(w+4)<<""<<"value   = U implies walk segment length value of 2^((Lodd + 1)/2 - laevus - 1)"<<endl;
    cout<<setw(w+8)<<""<<"(2 GB memory is required when (Lodd + 1)/2 - laevus  = 53)"<<endl;
    cout<<setw(w)<<"-seed s1,s2,s3";
    cout<<"... three **short** integers separated by commas, to initialize"<<endl;
    cout<<setw(w+4)<<""<<"the random number generator; default = ";
    cout<<"1,2,3."<<endl;
    cout<<setw(w)<<"-coordInit bitString";
    cout<<"... initial binary string,"<<endl;
    cout<<setw(w+4)<<""<<"default = a random binary string, determined by seed."<<endl;
    cout<<setw(w)<<"-v";
    cout<<"... turn on the verbose output."<<endl;
    #ifndef NDEBUG
    cout<<setw(w)<<"-trace";
    cout<<"... turn on the trace output."<<endl;
    cout<<setw(w)<<"-walk";
    cout<<"... turn on the walk output."<<endl;
    #endif
    cout<<setw(w)<<"-help";
    cout<<"... display this help and exit."<<endl;
    cout<<"------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Copyright 2014"<<endl;
    cout<<"*  Borko Bošković, Franc Brglez, and Janez Brest"<<endl;
    cout<<"*  Paper: Low-Autocorrelation Binary Sequences: On Improved"<<endl;
    cout<<"*    Merit Factors and Runtime Predictions to Achieve Them"<<endl;
    cout<<"*  Solver: Self-Avoiding Walk"<<endl;
    cout<<"*  arXiv:1406.5301 (http://arxiv.org/abs/1406.5301)"<<endl;
    cout<<"*  Version "<<VERSION<<endl;
    cout<<"*  Compile date "<<__DATE__<<endl;
    exit(0);
}

int main(int argc, char *argv[]){
    unsigned int uiValue;
    double dValue;
    char cValue;
    SAW saw;
    short unsigned int seed[3] = {1,2,3};
    saw.setRuntimeLmt(300);
    saw.setSeed(seed);
    saw.setWalkLength(8);
    saw.setLaevus(0);
    saw.setL(0);
    saw.setValueTarget(0);

    stringstream ss(stringstream::in | stringstream::out);
    string arg, args;
    for(int i = 1; i < argc; i++) args += string(argv[i]) + " ";
    saw.setArgs(args);

    if(args.find("-help") != string::npos || argc == 1){
        help(argv[0]);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        if(!ss.good()){
            cerr<<"Error: wrong value of parameter "<<argv[i-2]<<endl;
            exit(1);
        }
        arg = argv[i];
        if (i == 1) {
            ss<<argv[i]<<" ";
            ss>>uiValue;
            saw.setL(uiValue);
        }
        else if (arg == "-valueTarget") {
            ss<<argv[++i]<<" ";
            ss>>uiValue;
            saw.setValueTarget(uiValue);
        }
        else if (arg == "-runtimeLmt") {
            ss<<argv[++i]<<" ";
            ss>>uiValue;
            saw.setRuntimeLmt(uiValue);
        }
        else if(arg == "-seed") {
            ss<<argv[++i]<<" ";
            ss>>seed[0];
            ss>>cValue;
            ss>>seed[1];
            ss>>cValue;
            ss>>seed[2];
            saw.setSeed(seed);
        }
        else if (arg == "-walkSegmCoef") {
            ss<<argv[++i]<<" ";
            ss>>dValue;
            if(!ss.good()){
            	ss.ignore();
                ss.clear();
		char wl;
		ss>>wl;
		if(wl == 'U' || wl == 'u')
			saw.setWalkLength(0);
		else{
			cerr<<"Wrong walkLenght argument:"<<wl<<endl;	
			exit(1);
		}
            }
            else saw.setWalkLength(dValue);
        }
        else if (arg == "-laevus") {
            ss<<argv[++i]<<" ";
            ss>>uiValue;
            saw.setLaevus(uiValue);
        }
        else if (arg == "-v") {
            saw.setVerbose();
        }
        #ifndef NDEBUG
        else if (arg == "-trace") {
            saw.setTrace();
        }
        else if (arg == "-walk") {
            saw.setWalk();
        }
        #endif
        else if (arg == "-coordInit"){
            char coordInit[MAX_S_SIZE+1];
            ss<<argv[++i]<<" ";
            ss>>coordInit;
            saw.setCoordInit(coordInit);
        }
        else{
            cerr<<"Wrong parameter:"<<arg<<endl;
            exit(1);
        }
    }

    if(saw.getValueTarget() == 0 && saw.getL() < MAX_BEST_L)
        saw.setValueTarget(bestT[saw.getL()]);

    saw.run();
    return 0;
}
