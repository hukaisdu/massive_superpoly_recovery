#ifndef __GRAIN_H__
#define  __GRAIN_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

void BackExpandPolynomial( int rounds, vector<bitset<544>> & terms );

void SecondBackExpandPolynomial( int rounds, const bitset<288 + 256> & last,
		vector<bitset<288 + 256>> & terms, int threads );

STATUS MidSolutionCounter( int rounds, bitset<128> cube, const bitset<544> & last, map<bitset<128>, 
		int, CMPS<128>> & counterMap,  float time, int threads );

#endif
