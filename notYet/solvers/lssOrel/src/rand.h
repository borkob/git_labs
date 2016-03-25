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

#ifndef RAND_H
#define RAND_H

#include <iostream>
#include <cstdlib>
#include <cstdint>
using namespace std;

class Rand {
public:
    Rand (const unsigned short seed[3]) {
        for (unsigned int i=0; i<3; i++ ) this->seed[i] = seed[i];
    }
    inline int int32(const int n) { return (nrand48(seed)%n); }
    inline uint64_t int64() {
        return uint64_t(nrand48(seed))
                | (uint64_t(nrand48(seed))<<31)
                | (uint64_t(nrand48(seed))<<62);
    }
    inline double doubleValue(const double lbound, const double ubound ) {
        return lbound + (ubound - lbound) * erand48(seed);
    }
private:
    unsigned short seed[3];
};

#endif // RAND_H
