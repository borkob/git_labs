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

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "sequence.h"
#include "rand.h"
#include <stdint.h>
#include <cstring>

#ifndef MAX_L
#define MAX_L 320
#endif

#if MAX_L%64
#error MAX_L%64 != 0!
#endif

#define MAX_S_SIZE (MAX_L/2)

class Sequence
{
public:
    Sequence();
    static void setL(const unsigned int L, Rand * rand);
    static void setLaevus(const unsigned int laevus);
    #ifndef NDEBUG
    unsigned int evaluate();
    #endif
    void random(Rand * rand);
    void random(Rand * rand, const unsigned int n);
    void random(Rand * rand, const double biasW);
    void setS(char S[MAX_S_SIZE+1], Rand * rand);
    friend ostream& operator<<(ostream& out, const Sequence individual);
    friend istream& operator>>(istream& in, Sequence& individual);
    void neighborhoodStructure();
    uint64_t updateNeighborhoodStructure(const unsigned int bit, const unsigned int eval);
    unsigned int flipEval(const unsigned int bit) const;
    #ifndef NDEBUG
    void check();
    #endif
    string toString() const ;
    string toSkewString() const ;

    inline unsigned int getEval() const { return eval; }
    inline void setEval(const unsigned int eval)  { this->eval = eval; }
    inline uint64_t getHKey() const { return hKey; }
    inline static uint64_t getCntProbe() { return cntProbe; }
    inline static void resetCntProbe() { cntProbe = 0; }
    inline static void addCntProbe(unsigned int v){ cntProbe += v; }

    inline int getS(const unsigned int i) const{
        #ifndef NDEBUG
        if(i >= D){
            cerr<<"Wrong index (getS):"<<i<<endl;
            exit(1);
        }
        #endif
        return s[i];
    }
    inline void setS(const unsigned int i, const int value) {
        #ifndef NDEBUG
        if(i >= D){
            cerr<<"Wrong index (setS):"<<i<<endl;
            exit(1);
        }
        #endif
        s[i] = value;
    }
    #ifndef NDEBUG
    inline void flip(const unsigned int i){
        if(i >= D){
            cerr<<"Wrong index (flip):"<<i<<endl;
            exit(1);
        }
        hKey ^= rkey[i];
        s[i] = -s[i];
    }
    #endif
    inline void flip(const unsigned int i, const unsigned int eval){
        #ifndef NDEBUG
        if(i >= D){
            cerr<<"Wrong index (flip):"<<i<<endl;
            exit(1);
        }
        #endif
        hKey ^= rkey[i];
        s[i] = -s[i];
        this->eval = eval;
    }
    inline uint64_t getFlipedKey(const unsigned int i){
        #ifndef NDEBUG
        if(i >= D){
            cerr<<"Wrong index (flip):"<<i<<endl;
            exit(1);
        }
        #endif
        return hKey ^ rkey[i];
    }

private:
    int s[MAX_S_SIZE];
    unsigned int eval;
    uint64_t hKey;

    static uint64_t rkey[MAX_S_SIZE];
    static unsigned int L;
    static unsigned int L2;
    static unsigned int D;
    static unsigned int laevus;
    static uint64_t cntProbe;
    static int Tau[MAX_S_SIZE][MAX_S_SIZE];
    static int C[MAX_S_SIZE];
};

#endif // SEQUENCE_H
