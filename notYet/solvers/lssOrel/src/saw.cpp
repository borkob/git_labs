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

#include <sys/utsname.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <cmath>

const unsigned int SAW::w = 22;

SAW::SAW(){
    runtimeLmt = D = laevus = valueTarget = L = 0;
    walkLength = 0;
    walkLengthFactor = 8;
    rand = NULL;
    seed[0] = seed[1] = seed[2] = 0;
    verbose = false;
    coordInit[0] = '\0';

    utsname un;
    uname(&un);
    systemInfo = un.machine;
    systemInfo += " ";
    systemInfo += un.release;
    systemInfo += " ";
    systemInfo += un.sysname;
    progName = NAME;

    #ifndef NDEBUG
    trace = false;
    #endif
}

bool SAW::stoppingCriterion(){
    if(runtimeLmt != 0){
        cuurentTime = double(clock() - startTime) / CLOCKS_PER_SEC;
        if(cuurentTime >= (double)runtimeLmt) return true;
    }
    if(best.getEval() <= valueTarget) return true;
    return false;
}

void SAW::start(){
    if(L >= MAX_L || L == 0 || L%2 == 0){
        cerr<<"Error: oddL = "<<L<<"!"<<endl;
        exit(1);
    }
    if(laevus >= D){
        cerr<<"Error: laevus = "<<laevus<<"!"<<endl;
        exit(1);
    }

    if(walkLengthFactor == 0) walkLength = (unsigned long long) pow(2,D-laevus-1);
    else walkLength = walkLengthFactor*(D-laevus);

    rand = new Rand(seed);
    Sequence::resetCntProbe();
    Sequence::setL(L,rand);
    Sequence::setLaevus(laevus);
    startTime = clock() - 0.001;
    restart = 0;
    walkLengthTot = 0;
    //path.reserve(walkLength<<4);
    #ifndef NDEBUG
    step = 0;
    #endif
}

void SAW::stop(){
    printInfo();
    printResults();
    if(rand) delete rand;
    #ifndef NDEBUG
    if(walk) fileStream.close();
    #endif
}

void SAW::printHeader(ostream& stream, const bool outputFile){
    time_t now = time(0);
    if(outputFile) now = time(0);
    ctime(&now);
    stream<<left;
    const string solverName = progName.substr(0,progName.find("_"));
    stringstream ss(stringstream::in | stringstream::out);
    if(L < 10) ss<<"n=00";
    else if(L < 100) ss<<"n=0";
    else ss<<"n=";
    ss<<L<<"-B-"<<laevus<<"_"<<solverName<<"_"<<walkLength<<"-"<<trial<<"-walk0.txt";
    const string outFileName = ss.str();
    stream<<setw(w)<<"# "<<outFileName<<endl;
    stream<<"# -->>>------------------------------------------------------------";
    stream<<"------"<<endl;
    stream<<setw(w)<<"# date"<<ctime(&now);
    stream<<setw(w)<<"# functionName"<<"labs"<<endl;
    stream<<setw(w)<<"# instanceDef"<<L<<endl;
    stream<<setw(w)<<"# laevus"<<laevus<<endl;
    stream<<setw(w)<<"# isSkew"<<1<<endl;
    stream<<setw(w)<<"# solverName"<<solverName<<endl;
    stream<<setw(w)<<"# coordinateType"<<"B"<<endl;
    stream<<setw(w)<<"# nDim"<<D-laevus<<endl;
    if(valueTarget>0){
        stream<<setw(w)<<"# valueTarget"<<valueTarget<<endl;
        stream<<setw(w)<<"# meritTarget"<<setprecision(4)<<fixed<<(L*L)/(2.0*valueTarget)<<endl;
    }
    else{
        stream<<setw(w)<<"# valueTarget"<<"NA"<<endl;
        stream<<setw(w)<<"# meritTarget"<<"NA"<<endl;
    }
    stream<<setw(w)<<"# progName"<<progName<<endl;
    stream<<setw(w)<<"# progVersion"<<VERSION<<endl;
    stream<<setw(w)<<"# compile date"<<__DATE__<<endl;
    stream<<setw(w)<<"# args"<<args<<endl;
    stream<<setw(w)<<"# system"<<systemInfo<<endl;
    stream<<setw(w)<<"# algorithm"<<"Self-avoiding walk"<<endl;
    if(runtimeLmt != 0) stream<<setw(w)<<"# runtimeLmt"<<runtimeLmt<<endl;
    else stream<<setw(w)<<"# runtimeLmt"<<"infinite"<<endl;
    stream<<setw(w)<<"# walkSegmLmt"<<walkLength<<endl;
    stream<<setw(w)<<"# seed"<<seed[0]<<","<<seed[1]<<","<<seed[2]<<endl;
    if(coordInit[0] == '\0') stream<<setw(w)<<"# coordInit"<<"NA"<<endl;
    else stream<<setw(w)<<"# coordInit"<<coordInit<<endl;
    stream<<"# -->>>------------------------------------------------------------";
    stream<<"------"<<endl;
    #ifndef NDEBUG
    if(walk){
        stream<<"step\tchain\tcoord\tvalue\tcomment"<<endl;
        if(outputFile) return;
        fileStream.open(outFileName.c_str());
        if(!fileStream.is_open()){
            cerr<<"Error: file:"<<outFileName<<endl;
            exit(1);
        }
        printHeader(fileStream,true);
    }
    #endif
    if(verbose){
        stream<<setw(15)<<"# cntProbe";
        stream<<setw(13)<<"cntRestart";
        stream<<setw(7)<<"value";
        stream<<setw(8)<<"merit";
        stream<<setw(11)<<"runtime";
        stream<<setw(13)<<"speed";
        stream<<setw(20)<<"coord"<<endl;
    }
}

void SAW::printInfo(){
    if(!verbose) return;
    cuurentTime = double(clock() - startTime) / CLOCKS_PER_SEC;
    cout<<"# "<<setw(13)<<scientific<<(double)Sequence::getCntProbe();
    cout<<setw(13)<<scientific<<(double)restart;
    cout<<setw(7)<<best.getEval();
    cout<<setw(8)<<setprecision(4)<<fixed<<(L*L)/(2.0*best.getEval());
    cout<<setw(11)<<setprecision(2)<<fixed<<cuurentTime<<setprecision(6);
    cout<<setw(13)<<scientific<<Sequence::getCntProbe()/cuurentTime;
    cout<<best<<endl;
}

void SAW::printResults() {
    #ifndef NDEBUG
    if(trace){
        cout<<"walkLength=0 cntRestart="<<restart-1<<" cntProbe="<<Sequence::getCntProbe();
        cout<<" pivot="<<best<<":"<<best.getEval()<<"\t# PIVOT (TARGET) **END**"<<endl;
    }
    if(walk){
        if(restart == 0) restart ++;
        cout<<step<<"\t"<<restart<<"\t"<<best<<"\t"<<best.getEval();
        cout<<"\t# PIVOT (TARGET) **END**"<<endl;
        fileStream<<step<<"\t"<<restart<<"\t"<<best<<"\t"<<best.getEval();
        fileStream<<"\t# PIVOT (TARGET) **END**"<<endl;
        step ++;
        return;
    }
    #endif

    cuurentTime = double(clock() - startTime) / CLOCKS_PER_SEC;
    Sequence init;
    init.setS(coordInit, NULL);
    time_t now = time(0);
    ctime(&now);

    cout<<left;
    cout<<"# -->>>------------------------------------------------------------";
    cout<<"------"<<endl;
    cout<<setw(w)<<"# date"<<ctime(&now);
    cout<<setw(w)<<"# args"<<args<<endl;
    cout<<setw(w)<<"# version"<<VERSION<<endl;
    cout<<setw(w)<<"# compile date"<<__DATE__<<endl;
    cout<<setw(w)<<"# system"<<systemInfo<<endl;
    cout<<setw(w)<<"# algorithm"<<"Self-avoiding walk"<<endl;
    cout<<setw(w)<<"instanceDef"<<L<<endl;
    cout<<setw(w)<<"nDim"<<D-laevus<<endl;
    cout<<setw(w)<<"progName"<<progName<<endl;
    cout<<setw(w)<<"progVersion"<<VERSION<<endl;
    if(valueTarget>0){
        cout<<setw(w)<<setprecision(4)<<fixed<<"meritTarget"<<(L*L)/(2.0*valueTarget)<<endl;
        cout<<setw(w)<<"valueTarget"<<valueTarget<<endl;
    }
    else{
        cout<<setw(w)<<"meritTarget"<<"NA"<<endl;
        cout<<setw(w)<<"valueTarget"<<"NA"<<endl;
    }
    cout<<setw(w)<<"valueBest"<<best.getEval()<<endl;
    cout<<setw(w)<<"meritBest"<<(L*L)/(2.0*best.getEval())<<endl;
    if(best.getEval() > valueTarget) cout<<setw(w)<<"targetReached"<<0<<endl;
    else if(best.getEval() == valueTarget) cout<<setw(w)<<"targetReached"<<1<<endl;
    else cout<<setw(w)<<"targetReached"<<2<<endl;
    cout<<setw(w)<<"isCensored"<<(runtimeLmt != 0 && runtimeLmt<=cuurentTime)<<endl;
    if(runtimeLmt != 0) cout<<setw(w)<<"runtimeLmt"<<runtimeLmt<<endl;
    else cout<<setw(w)<<"runtimeLmt"<<"NA"<<endl;
    cout<<setw(w)<<"runtime"<<setprecision(2)<<fixed<<cuurentTime<<setprecision(6)<<endl;
    cout<<setw(w)<<"cntProbe"<<Sequence::getCntProbe()<<endl;
    cout<<setw(w)<<"cntRestart"<<restart<<endl;
    cout<<setw(w)<<"walkLength"<<walkLengthTot<<endl;
    cout<<setw(w)<<"speed"<<(uint64_t)(Sequence::getCntProbe()/cuurentTime)<<endl;
    cout<<setw(w)<<"functionParameters"<<laevus<<endl;
    if(walkLengthFactor != 0)
	    cout<<setw(w)<<"walkSegmCoef"<<setprecision(2)<<walkLengthFactor<<setprecision(6)<<endl;
    else
	    cout<<setw(w)<<"walkSegmCoef "<<"NA"<<endl;
    cout<<setw(w)<<"seedFirst"<<seed[0]<<","<<seed[1]<<","<<seed[2]<<endl;
    cout<<setw(w)<<"coordInit"<<init.toSkewString()<<endl;
    cout<<setw(w)<<"coordInitFull"<<init.toString()<<endl;
    cout<<setw(w)<<"coordBest"<<best.toSkewString()<<endl;
    cout<<setw(w)<<"coordBestFull"<<best.toString()<<endl;
}

void SAW::run(){
    start();
    trial.setS(coordInit,rand);
    trial.neighborhoodStructure();
    best = trial;
    #ifndef NDEBUG
    if(verbose || walk) printHeader(cout);
    #else
    if(verbose){
        printHeader(cout);
        printInfo();
    }
    #endif
    while(!stoppingCriterion()){
        saw();
        if(trial.getEval() < best.getEval()){
            best = trial;
            printInfo();
            if(best.getEval() <= valueTarget) break;
        }
        restart++;
	walkLengthTot++;
        trial.random(rand);
        trial.neighborhoodStructure();
    }
    #ifndef NDEBUG
    if(walk && Sequence::getCntProbe() == 1){
        cout<<step<<"\t"<<restart<<"\t"<<trial<<"\t"<<trial.getEval();
        fileStream<<step<<"\t"<<restart<<"\t"<<trial<<"\t"<<trial.getEval();
        if(trial.getEval() <= valueTarget){
            cout<<"\t# targetReached 1 cntProbe "<<Sequence::getCntProbe()<<endl;
            fileStream<<"\t# targetReached 1 cntProbe "<<Sequence::getCntProbe()<<endl;
        }
        else{
            cout<<"\t#"<<endl;
            fileStream<<"\t#"<<endl;
        }
        step ++;
    }
    #endif
    stop();
}

void SAW::saw() {
    unsigned int bestMoveEval, fCurrent, fLocalBest, i, iBefore;
    unsigned int bestMovesBits[D], numBestMoves;
    path.clear();
    Sequence current = trial;
    fLocalBest = current.getEval();
    #ifndef NDEBUG
    int tmpProbe = 0;
    current.check();
    if(trace){
        cout<<"walkLength=0 cntRestart="<<restart<<" cntProbe=";
        cout<<Sequence::getCntProbe()<<" pivotInit="<<current<<":";
        cout<<current.getEval()<<endl;
    }
    if(walk){
        cout<<step<<"\t"<<restart<<"\t"<<current<<"\t"<<current.getEval();
        fileStream<<step<<"\t"<<restart<<"\t"<<current<<"\t"<<current.getEval();
        if(current.getEval() <= valueTarget){
            cout<<"\t# targetReached 1 cntProbe "<<Sequence::getCntProbe()<<endl;
            fileStream<<"\t# targetReached 1 cntProbe "<<Sequence::getCntProbe()<<endl;
        }
        else{
            cout<<"\t#"<<endl;
            fileStream<<"\t#"<<endl;
        }
        step ++;
    }
    #endif
    iBefore = numeric_limits<unsigned int>::max();
    for(unsigned long long int k=walkLength; k>0; k--){
        path.insert(current.getHKey());
        bestMoveEval = numeric_limits<unsigned int>::max();
        numBestMoves = 0;
        for(i=laevus; i<D; i++){
            if(i == iBefore) continue;
            fCurrent = current.flipEval(i);;
            #ifndef NDEBUG
            if(trace){
                tmpProbe++;
                current.flip(i);
                cout<<"\tprobe="<<Sequence::getCntProbe()+tmpProbe<<" flipped bit = ";
                cout<<(i+1-laevus)<<", "<<current<<":"<<fCurrent<<endl;
                current.flip(i);
            }
            #endif
            if(fCurrent < bestMoveEval){
                if(path.count(current.getFlipedKey(i)) > 0) continue;
                bestMovesBits[0] = i;
                numBestMoves = 1;
                bestMoveEval = fCurrent;
            }
            else if(fCurrent == bestMoveEval){
                if(path.count(current.getFlipedKey(i)) > 0) continue;
                bestMovesBits[numBestMoves] = i;
                numBestMoves++;
            }
        }
        #ifndef NDEBUG
        tmpProbe = 0;
        #endif
        if(k != walkLength) Sequence::addCntProbe(D-laevus-1);
        else Sequence::addCntProbe(D-laevus);
        if(numBestMoves > 0){
            walkLengthTot++;
            if(numBestMoves == 1) iBefore = bestMovesBits[0];
            else iBefore = bestMovesBits[rand->int32(numBestMoves)];
            current.updateNeighborhoodStructure(iBefore,bestMoveEval);
            #ifndef NDEBUG
            current.check();
            if(trace){
                cout<<"walkLength="<<(walkLength-k+1)<<" cntRestart="<<restart<<" cntProbe=";
                cout<<Sequence::getCntProbe()<<" bit="<<(iBefore+1)<<" pivot=";
                cout<<current<<":"<<current.getEval()<<" bestMoves "<<numBestMoves<<endl;
            }
            if(walk){
                cout<<step<<"\t"<<restart<<"\t"<<current<<"\t"<<current.getEval();
                fileStream<<step<<"\t"<<restart<<"\t"<<current<<"\t"<<current.getEval();
                if(current.getEval() <= valueTarget){
                    cout<<"\t# targetReached 1 cntProbe "<<Sequence::getCntProbe()<<endl;
                    fileStream<<"\t# targetReached 1 cntProbe "<<Sequence::getCntProbe()<<endl;
                }
                else{
                    cout<<"\t#"<<endl;
                    fileStream<<"\t#"<<endl;
                }
                step ++;
            }
            #endif
            if(bestMoveEval < fLocalBest) {
                fLocalBest = bestMoveEval;
                trial = current;
                if(fLocalBest <= valueTarget) return;
            }
        }
        else return;
    }
}
