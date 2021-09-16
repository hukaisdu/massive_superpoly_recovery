#ifndef __GRAIN_H__
#define  __GRAIN_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

int BackExpandPolynomial( int rounds, vector<bitset<256>> & term );

int SecondBackExpandPolynomial( int rounds, const bitset<256> & last,
		vector<bitset<256> > & term, int threads  );
STATUS MidSolutionCounter( int rounds, const bitset<96> & cube, const
		bitset<256> & last, map<bitset<128>, int, CMPS<128>> & counterMap, float
		time, int thread );

#endif
