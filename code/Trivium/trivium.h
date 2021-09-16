#ifndef __TRIVIUM_H__
#define __TRIVIUM_H__
#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"

void triviumCore( GRBModel&, vector<GRBVar>&, int, int, int, int, int );

int SecondBackExpandPolynomial( int, const bitset<288> &, vector<bitset<288> > &, int );

STATUS MidSolutionCounter( int, const bitset<80> &, const bitset<288> &,
		map<bitset<80>, int, CMPS<80>> &, float, int );

int BackExpandPolynomial( int, vector<bitset<288> > & );

#endif

