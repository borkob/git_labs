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

#ifndef SAW_H
#define SAW_H

#include "rand.h"
#include "sequence.h"

#include <unordered_set>
#include <cstring>
#include <fstream>

class SAW
{
public:
    SAW();
    void run();

    inline void setL(const unsigned int L){
        if(L >= MAX_L){
            cerr<<"Error: L = "<<L<<"!"<<endl;
            exit(1);
        }
        this->L = L;
        D=(L+1)/2;
    }

    inline void setWalkLength(const double walkLength) {
        walkLengthFactor = walkLength;
    }
    inline void setLaevus(const unsigned int laevus){ this->laevus = laevus; }
    inline void setName(const string & name){ this->progName = name; }
    inline void setArgs(const string & args){ this->args = args; }
    inline void setValueTarget(const unsigned int targetValue) { this->valueTarget = targetValue; }
    inline void setRuntimeLmt(const unsigned int runtimeLmt) { this->runtimeLmt = runtimeLmt; }
    inline void setSeed(const unsigned short seed[3]) { memcpy(this->seed,seed,3*sizeof(short)); }
    inline void setVerbose(){ verbose = true; }
    inline void setCoordInit(const char coordInit[MAX_S_SIZE+1]) { memcpy(this->coordInit,coordInit,MAX_S_SIZE+1); }

    inline unsigned int getL() const { return L; }
    inline unsigned int getLaevus() const { return laevus; }
    inline unsigned int getValueTarget() const { return valueTarget; }
    inline unsigned int getRuntimeLmt() const { return runtimeLmt; }
    inline unsigned short int getSeed(const unsigned int i) const { return seed[i]; }
    inline unsigned int getD() const { return D; }
    #ifndef NDEBUG
    inline void setTrace() { trace = true;};
    inline void setWalk() { walk = true;};
    #endif
private:
    bool stoppingCriterion();
    void start();
    void stop();
    void printHeader(ostream& stream, const bool outputFile=false);
    void printInfo();
    void printResults();
    void saw();

    Sequence trial;
    Sequence best;

    unsigned int runtimeLmt;
    clock_t startTime;
    double cuurentTime;
    unsigned int restart;

    unsigned long long int walkLength;
    unsigned long long int walkLengthTot;
    double walkLengthFactor;
    unsigned int L;
    unsigned int D;
    unsigned int laevus;
    unsigned int valueTarget;
    unsigned short seed[3];

    static const unsigned int w;
    #ifndef NDEBUG
    bool trace;
    bool walk;
    unsigned int step;
    ofstream fileStream;
    #endif

    char coordInit[MAX_S_SIZE+1];

    string progName;
    string args;
    string systemInfo;

    Rand * rand;

    bool verbose;
    unordered_set<uint64_t> path;

};

#endif // TABUSEARCH_H
