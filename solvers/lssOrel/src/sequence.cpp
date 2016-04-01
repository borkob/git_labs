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

#include "sequence.h"

#include <limits>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

unsigned int Sequence::L;
unsigned int Sequence::L2;
unsigned int Sequence::D;
unsigned int Sequence::laevus;
int Sequence::Tau[MAX_S_SIZE][MAX_S_SIZE];
int Sequence::C[MAX_S_SIZE];
uint64_t Sequence::cntProbe = 0;
uint64_t Sequence::rkey[MAX_S_SIZE];

Sequence::Sequence(){
    memset(s,0,MAX_S_SIZE*sizeof(int));
    hKey = 0;
    eval =  numeric_limits<unsigned int>::max();
}

void Sequence::setL(const unsigned int L, Rand * rand){
    if(L >= MAX_L || L == 0 || L%2 == 0){
        cerr<<"Error: L = "<<L<<"!"<<endl;
        exit(1);
    }
    D = (L+1)/2;
    Sequence::L = L;
    L2 = L/2;
    for(unsigned int i=0; i<D; i++) rkey[i] = rand->int64();
}

void Sequence::setLaevus(const unsigned int laevus){
    if(laevus >= MAX_S_SIZE){
        cerr<<"Error: z = "<<laevus<<"!"<<endl;
        exit(1);
    }
    Sequence::laevus = laevus;
}

#ifndef NDEBUG
unsigned int Sequence::evaluate() {
    int ck;
    unsigned int n = D-1;
    int ss[L];
    memcpy(ss,s,D*sizeof(int));
    for(unsigned int i=1; i<=n; i++){
        if(i%2==0) ss[n+i] = ss[n-i];
        else ss[n+i] = -ss[n-i];
    }
    eval = 0;
    for (unsigned int k=2; k<=L-1; k+=2) {
        ck=0;
        for (unsigned int i=0; i<=L-k-1; i++) ck += ss[i]*ss[i+k];
        eval+=ck*ck;
    }
    cntProbe++;
    return eval;
}
#endif

void Sequence::random(Rand * rand){
    hKey = 0;
    for(unsigned int i=laevus; i<D; i++){
        if(rand->doubleValue(0,1)<0.5) s[i] = 1; else s[i] = -1;
        if(s[i] == 1) hKey ^= rkey[i];
    }
}

void Sequence::random(Rand * rand, const unsigned int n){
    static const unsigned int nDim = D - laevus;
    unsigned int index, changes[n];
    bool found;
    for(unsigned int i=0; i<n; i++){
        while(true){
            index = laevus + rand->int32(nDim);
            found = false;
            for(unsigned int j=0; j<i; j++){
                if(changes[j] == index){
                    found = true;
                    break;
                }
            }
            if(!found){
                changes[i] = index;
                break;
            }
        }
        s[index] *= -1;
    }
}

void Sequence::setS(char S[MAX_S_SIZE+1], Rand * rand){
    if(S[0] != '\0' && (D-laevus) != strlen(S)){
        cerr<<"Wrong char sequence: "<<S;
        cerr<<" - z("<<laevus<<") + S("<<strlen(S)<<") != Ls("<<D<<")"<<endl;
        exit(1);
    }
    if(S[0] == '\0'){
        for(unsigned int i=laevus; i<D; i++){
            if(rand->doubleValue(0,1)<0.5) S[i-laevus] = '1';
            else S[i-laevus] = '0';
        }
        S[D-laevus] = '\0';
    }
    hKey = 0;
    for(unsigned int i=0; i<D; i++){
        if(i < laevus) s[i] = -1;
        else if(S[i-laevus]=='1') s[i] = 1;
        else if(S[i-laevus]=='0') s[i] = -1;
        else{
            cerr<<"Error: wrong char sequence:"<<S[i-laevus]<<endl;
            exit(1);
        }
        if(s[i] == 1) hKey ^= rkey[i];
    }
}

string Sequence::toString() const {
    string str = "";
    unsigned int n = D-1;
    int ss[L];
    memcpy(ss,s,D*sizeof(int));
    for(unsigned int i=1; i<=n; i++){
        if(i%2==0) ss[n+i] = ss[n-i];
        else ss[n+i] = -ss[n-i];
    }

    for(unsigned int i=0; i<L; i++){
        if(ss[i] == 1) str+="1";
        else str+="0";
    }
    return str;
}

string Sequence::toSkewString() const {
    string str = "";
    for(unsigned int i=0; i<D; i++){
        if(s[i] == 1) str+="1";
        else str+="0";
    }
    return str;
}

ostream& operator<<(ostream& out, Sequence seq) {
    for(unsigned int i=Sequence::laevus; i<Sequence::D; i++){
        if(seq.s[i] == 1) out<<"1";
        else if(seq.s[i] == -1) out<<"0";
        else{
            cerr<<"wrong binary string:"<<seq.s[i]<<endl;
            exit(1);
        }
    }
    return out;
}

istream& operator>>(istream& in, Sequence& seq) {
    char x;
    for(unsigned int i=0; i< Sequence::laevus; i++)
        seq.s[i] = 0;
    for(unsigned int i=Sequence::laevus; i<seq.D; i++){
        in>>x;
        if(x == '1') seq.s[i] = 1;
        else seq.s[i] = -1;
    }
    return in;
}

/**
 * In order to make the performance comparisons of solvers 
 * lssOrel, lssMAta, and lssRRts as described in our paper
 *   Borko Bošković, Franc Brglez, and Janez Brest
 *   Low-Autocorrelation Binary Sequences: on the Performance
 *   of Memetic-Tabu and Self-Avoiding Walk Solvers
 *   http://arxiv.org/abs/1406.5301
 *
 * as accurate as possible, our solver lssOrel uses the same methods
 * - neighborhoodStructure
 * - updateNeighborhoodStructure
 * - flipEval
 *
 * that are being used in solvers lssMAta, and lssRRts for the incremental 
 * value calculation. These two solvers are the instrumented version of 
 * the solver MAts_{skew} described in the paper
 *   José E. Gallardo, Carlos Cotta, and Antonio J. Fernández.
 *   Finding low autocorrelation binary sequences with memetic algorithms.
 *   Appl. Soft Comput., 9(4):1252–1262, September 2009.
 *
 * and made available in our research by the authors.

 * With this approach, we can accurately compare merits of concepts 
 * such as self-avoiding walk, tabu search, and the evolutionary component. 
 * Again, all solvers, lssOrel, lssMAta, and lssRRts use the same code for 
 * the incremental value calculation. Within lssOrel, these calculations are 
 * implemented inside the ‘Sequence class’, and are invoked within the 
 * self-avoiding walk algorithm (SAW class).
*/

#define skew_s(i) ((i) < L/2 ? s[i] : (((L-1)-(i))%2==((L/2)%2) ? s[(L-1)-(i)] : -s[(L-1)-(i)] ) )

void Sequence::neighborhoodStructure() {
    for(unsigned int i=0; i<L2; i++){
        for(unsigned int j=0; j<L2-i; j++){
            Tau[i][j] = 4*s[j]*skew_s(j+(2*(i+1)));
        }
    }
    eval = 0;
    for(unsigned int i=0; i<L2; i++) {
        C[i] = Tau[i][L2-i-1]/4;
        for(unsigned int j=0; j<L2-i-1; j++)
            C[i] += Tau[i][j]/2;
        eval += C[i]*C[i];
    }
    cntProbe++;
}

uint64_t Sequence::updateNeighborhoodStructure(const unsigned int bit, const unsigned int eval){
    int i, val, *ptrTau, *ptrC, *ptrIEnd;
    for(ptrC=&C[0], ptrTau=&Tau[0][bit-2], ptrIEnd=ptrC+bit/2;
        ptrC<ptrIEnd;
        ptrTau+=(MAX_S_SIZE-2),ptrC++) {
            val = *ptrTau;
            *ptrTau = -val;
            *ptrC -= val;
    }
    if(L-1-bit != bit) {
        for(i=L2-bit, ptrC=&C[i], ptrIEnd=ptrC+(bit/2),ptrTau=&Tau[i][bit-2];
            ptrC<ptrIEnd;
            ptrTau+=(MAX_S_SIZE-2),ptrC++) {
                val = *ptrTau;
                *ptrTau = -val;
                *ptrC -= val;
        }
    }
    for(ptrC=&C[0], ptrIEnd=ptrC+L2-(bit+1), ptrTau=&Tau[0][bit];
        ptrC<ptrIEnd;
        ptrTau+=MAX_S_SIZE,ptrC++) {
            val = *ptrTau;
            *ptrTau = -val;
            *ptrC -= val;
    }
    hKey ^= rkey[bit];
    s[bit] = -s[bit];
    this->eval = eval;
    return hKey;
}

unsigned int Sequence::flipEval(const unsigned int bit) const{
    static const unsigned int sizeC = D*sizeof(int);
    int i, *ptrTau, *ptrC, *ptrIEnd;
    int CC[D];
    memcpy(CC,C,sizeC);
    for(ptrC=&CC[0], ptrTau=&Tau[0][bit-2], ptrIEnd=ptrC+bit/2;
        ptrC<ptrIEnd;
        ptrTau+=(MAX_S_SIZE-2),ptrC++)
            *ptrC -= *ptrTau;
    if(L-1-bit != bit) {
        for(i=(L2)-bit, ptrC=&CC[i], ptrIEnd=ptrC+(bit/2),ptrTau=&Tau[i][bit-2];
            ptrC<ptrIEnd;
            ptrTau+=(MAX_S_SIZE-2),ptrC++)
                *ptrC -= *ptrTau;
    }
    for(ptrC=&CC[0], ptrIEnd=ptrC+(L2)-(bit+1), ptrTau=&Tau[0][bit];
        ptrC<ptrIEnd;
        ptrTau+=MAX_S_SIZE,ptrC++)
            *ptrC -= *ptrTau;
    int eval = 0;
    for(ptrC=&CC[0],ptrIEnd=ptrC+(L2); ptrC<ptrIEnd; ptrC++)
        eval += (*ptrC)*(*ptrC);
    return eval;
}

#ifndef NDEBUG
void Sequence::check() {
    uint64_t tmpHKey = 0;
    for(unsigned int i = laevus; i<D; i++){
        if(s[i] == 1) tmpHKey ^= rkey[i];
    }
    if(hKey != tmpHKey){
        cerr<<"Wrong hash key!"<<endl;
        exit(1);
    }
    unsigned int tmpEval = evaluate();
    cntProbe --;
    if(tmpEval != eval){
        cerr<<"Wrong eval!"<<endl;
        exit(1);
    }
    for(unsigned int i=0; i<laevus; i++){
        if(s[i] != -1){
            cerr<<"Wrong sequence (laevus="<<laevus<<"):"<<*this<<endl;
            exit(1);
        }
    }
}
#endif

void Sequence::random(Rand * rand, const double biasW){
    const unsigned int target = biasW * D;
    const unsigned int runsTarget = target + rand->int32(3) - 1;
    const unsigned int weightTarget = target + rand->int32(3) - 1;
    const unsigned int iterationMax = 999990;
    unsigned int iteration = 0, weight, runs;
    while(iteration < iterationMax){
        iteration ++;
        runs = weight = 0;
        for(unsigned int i=0; i<D; i++){
            s[i] = int(biasW + rand->doubleValue(0,1));
            if(s[i] == 1) weight++;
            else s[i] = -1;
            if(i > 0 && s[i] != s[i-1]) runs++;
        }
        if(runs == runsTarget && weight == weightTarget) break;
    }
}
