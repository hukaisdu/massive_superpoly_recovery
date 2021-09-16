#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<algorithm>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"

using namespace std;

const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

GRBVar LFSR( GRBModel & model, vector <GRBVar> & x )
{
	GRBVar a = model.addVar( 0, 1, 0, GRB_BINARY );
	GRBVar b = model.addVar( 0, 1, 0, GRB_BINARY );

	// copy
	model.addConstr( x[0] <= a + b );
	model.addConstr( x[0] >= a );
	model.addConstr( x[0] >= b );

	vector<GRBVar> t = x;
	for ( int i = 0; i < 127; i++ )
		x[i] = t[i+1];
	x[127] = a;
	return b;
}

void triviumCore(GRBModel& model, vector<GRBVar>& x, int i1, int i2, int i3, int i4, int i5)
{
    GRBVar y1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y2 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y3 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y4 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y5 = model.addVar(0, 1, 0, GRB_BINARY);

    GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);

    // z3 and z4 are not needed, since z3 = z4 = a
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);

    //copy
    model.addConstr(y1 <= x[i1]);
    model.addConstr(z1 <= x[i1]);
    model.addConstr(y1 + z1 >= x[i1]);

    //copy
    model.addConstr(y2 <= x[i2]);
    model.addConstr(z2 <= x[i2]);
    model.addConstr(y2 + z2 >= x[i2]);

    //copy
    model.addConstr(y3 <= x[i3]);
    model.addConstr(a <= x[i3]);
    model.addConstr(y3 + a >= x[i3]);
    
    //copy
    model.addConstr(y4 <= x[i4]);
    model.addConstr(a <= x[i4]);
    model.addConstr(y4 + a >= x[i4]);
    //XOR
    model.addConstr(y5 == x[i5] + a + z1 + z2);

    x[i1] = y1;
    x[i2] = y2;
    x[i3] = y3;
    x[i4] = y4;
    x[i5] = y5;
}

void SecondBackExpandPolynomial( int rounds, const bitset<288 + 256> & last,
		vector<bitset<288 + 256>> & terms, int threads )
{
    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, 2000000000 );
    //env.set(GRB_IntParam_PoolSolutions, 100 );
    env.set(GRB_IntParam_Threads, threads);

    // Create the model
    GRBModel model = GRBModel(env);


	// register
	vector<GRBVar> is(288);
	for ( int i = 0; i < 288; i++ )
		is[i] = model.addVar(0, 1, 0, GRB_BINARY );

	vector<GRBVar> ikey(128);
	for ( int i = 0; i < 128; i++ )
		ikey[i] = model.addVar(0, 1, 0, GRB_BINARY );

	vector<GRBVar> iiv(128);
	for ( int i = 0; i < 128; i++ )
		iiv[i] = model.addVar(0, 1, 0, GRB_BINARY );

	vector<GRBVar> s = is;
	vector<GRBVar> key = ikey;
	vector<GRBVar> iv = iiv;

	// update
	for ( int i = 0; i < rounds; i++ )
	{
		// LFSR
		GRBVar a = LFSR( model, key );
		GRBVar b = LFSR( model, iv );

        triviumCore(model, s, 65, 170, 90, 91, 92);
        triviumCore(model, s, 161, 263, 174, 175, 176);
        triviumCore(model, s, 242, 68, 285, 286, 287);

		GRBVar t1 = model.addVar(0, 1, 0, GRB_BINARY );
		GRBVar t3 = model.addVar(0, 1, 0, GRB_BINARY );

		model.addConstr( t1 == s[92] + b );
		model.addConstr( t3 == s[287] + a );

		s[92] = t1;
		s[287] = t3;

        vector<GRBVar> temp = s;
        for (int i = 0; i < 288; i++) 
            s[(i + 1) % 288] = temp[i];
	}

	// key + iv + 288
	for ( int i = 0; i < 128; i++ )
		model.addConstr( key[i] == last[i] );
	for ( int i = 0; i < 128; i++ )
		model.addConstr( iv[i] == last[i + 128] );
	for ( int i = 0; i < 288; i++ )
		model.addConstr( s[i] == last[i + 256] );

	model.optimize();
	map<bitset<256 + 288>, int, CMPS<256 + 288>> counterMap;
	
    if ( model.get(GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
		if ( solCount >= MAX )
		{
			cerr << "solCount value  is too big !" << endl;
			exit(-1);
		}

		bitset<544> start;
        for ( int i = 0; i < solCount; i++ )
        {
            model.set(GRB_IntParam_SolutionNumber, i );

            for ( int j = 0; j < 128; j++ ) 
                if ( round( ikey[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;

            for ( int j = 0; j < 128; j++ ) 
                if ( round( iiv[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j + 128] = 1;
                else 
                    start[j + 128] = 0;

            for ( int j = 0; j < 288; j++ ) 
                if ( round( is[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j + 256] = 1;
                else 
                    start[j + 256] = 0;

            counterMap[start]++;
        }
    }
	else
	{
		cout << "Infeasible" << endl;
		exit(-1);
	}
	for ( auto it : counterMap )
		if ( it.second % 2 == 1 )
			terms.push_back( it.first );
}

void BackExpandPolynomial( int rounds, vector<bitset<544>> & terms )
{
    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, 2000000000 );
    //env.set(GRB_IntParam_PoolSolutions, 100 );
    //env.set(GRB_IntParam_Threads, 2);

    // Create the model
    GRBModel model = GRBModel(env);

	// Create variables
	// register
	vector<GRBVar> is(288);
	for ( int i = 0; i < 288; i++ )
		is[i] = model.addVar(0, 1, 0, GRB_BINARY );

	vector<GRBVar> ikey(128);
	for ( int i = 0; i < 128; i++ )
		ikey[i] = model.addVar(0, 1, 0, GRB_BINARY );

	vector<GRBVar> iiv(128);
	for ( int i = 0; i < 128; i++ )
		iiv[i] = model.addVar(0, 1, 0, GRB_BINARY );

	vector<GRBVar> s = is;
	vector<GRBVar> key = ikey;
	vector<GRBVar> iv = iiv;


	// update
	for ( int i = 0; i < rounds; i++ )
	{
		// LFSR
		GRBVar a = LFSR( model, key );
		GRBVar b = LFSR( model, iv );

        triviumCore(model, s, 65, 170, 90, 91, 92);
        triviumCore(model, s, 161, 263, 174, 175, 176);
        triviumCore(model, s, 242, 68, 285, 286, 287);

		GRBVar t1 = model.addVar(0, 1, 0, GRB_BINARY );
		GRBVar t3 = model.addVar(0, 1, 0, GRB_BINARY );

		model.addConstr( t1 == s[92] + b );
		model.addConstr( t3 == s[287] + a );

		s[92] = t1;
		s[287] = t3;

        vector<GRBVar> temp = s;
        for (int i = 0; i < 288; i++) 
            s[(i + 1) % 288] = temp[i];
	}

    GRBLinExpr nk = 0;
    for ( int i = 0; i < 288; i++ )
        if ( (i == 65) || (i == 92) || (i == 161) || (i == 176) || (i == 242) || (i == 287))
            nk += s[i];
        else 
            model.addConstr( s[i] == 0 );
    model.addConstr( nk == 1 );

	for ( int i = 0; i < 128; i++ )
		model.addConstr( key[i] == 0 );

	for ( int i = 0; i < 128; i++ )
		model.addConstr( iv[i] == 0 );

	model.optimize();
	map<bitset<544>, int, CMPS<544>> counterMap;
	
    if ( model.get(GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
		if ( solCount >= MAX )
		{
			cerr << "solCount value  is too big !" << endl;
			exit(-1);
		}

		bitset<544> start;
        for ( int i = 0; i < solCount; i++ )
        {
            model.set(GRB_IntParam_SolutionNumber, i );

            for ( int j = 0; j < 128; j++ ) 
                if ( round( ikey[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;

            for ( int j = 0; j < 128; j++ ) 
                if ( round( iiv[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j + 128] = 1;
                else 
                    start[j + 128] = 0;

            for ( int j = 0; j < 288; j++ ) 
                if ( round( is[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j + 256] = 1;
                else 
                    start[j + 256] = 0;
            counterMap[start]++;
        }
    }
	else
	{
		cout << "Infeasible" << endl;
		exit(-1);
	}
	for ( auto it : counterMap )
		if ( it.second % 2 == 1 )
			terms.push_back( it.first );
}

STATUS MidSolutionCounter( int rounds, bitset<128> cube, const bitset<544> & last, map<bitset<128>, 
		int, CMPS<128>> & counterMap,  float time, int threads )
{
    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, 2000000000 );
    //env.set(GRB_IntParam_PoolSolutions, 100 );
    env.set(GRB_IntParam_Threads, threads);

    // Create the model
    GRBModel model = GRBModel(env);

	// Create variables
	vector<GRBVar> KEY(128);
	for ( int i = 0; i < 128; i++ )
		KEY[i] = model.addVar(0, 1, 0, GRB_BINARY );
	vector<GRBVar> IV(128);
	for ( int i = 0; i < 128; i++ )
		IV[i] = model.addVar(0, 1, 0, GRB_BINARY );
	
	// set cube
	for ( int i = 0; i < 128; i++ )
		if ( cube[i] == 1 )
			model.addConstr( IV[i] == 1 );
		else
			model.addConstr( IV[i] == 0 );

	// register
	vector<GRBVar> s(288);
	for ( int i = 0; i < 288; i++ )
		s[i] = model.addVar(0, 1, 0, GRB_BINARY );
	vector<GRBVar> key(128);
	for ( int i = 0; i < 128; i++ )
		key[i] = model.addVar(0, 1, 0, GRB_BINARY );
	vector<GRBVar> iv(128);
	for ( int i = 0; i < 128; i++ )
		iv[i] = model.addVar(0, 1, 0, GRB_BINARY );

	// distribute key
	for ( int i = 0; i < 128; i++ )
	{
		if ( i < 93 )
		{
			model.addConstr( KEY[i] <= s[i] + key[127 - i] );
			model.addConstr( KEY[i] >= s[i] );
			model.addConstr( KEY[i] >= key[127 - i] );
		}
		else
		{
			model.addConstr( KEY[i] == key[127 - i] );
		}
	}
	// distribute iv
	for ( int i = 0; i < 128; i++ )
	{
		model.addConstr( IV[i] <= iv[127 - i] + s[93 + i] );
		model.addConstr( IV[i] >= iv[127 - i] );
		model.addConstr( IV[i] >= s[93 + i] );
	}

	// constant
	model.addConstr( s[287] == 0 );

	// update
	for ( int i = 0; i < rounds; i++ )
	{
		// LFSR
		GRBVar a = LFSR( model, key );
		GRBVar b = LFSR( model, iv );

        triviumCore(model, s, 65, 170, 90, 91, 92);
        triviumCore(model, s, 161, 263, 174, 175, 176);
        triviumCore(model, s, 242, 68, 285, 286, 287);

		GRBVar t1 = model.addVar(0, 1, 0, GRB_BINARY );
		GRBVar t3 = model.addVar(0, 1, 0, GRB_BINARY );

		model.addConstr( t1 == s[92] + b );
		model.addConstr( t3 == s[287] + a );

		s[92] = t1;
		s[287] = t3;

        vector<GRBVar> temp = s;
        for (int i = 0; i < 288; i++) 
            s[(i + 1) % 288] = temp[i];
	}

	// key + iv + 288
	for ( int i = 0; i < 128; i++ )
		model.addConstr( key[i] == last[i] );
	for ( int i = 0; i < 128; i++ )
		model.addConstr( iv[i] == last[i + 128] );
	for ( int i = 0; i < 288; i++ )
		model.addConstr( s[i] == last[i + 256] );

    GRBLinExpr nv = 0;
    for ( int i = 0; i < 128; i++ )
		nv += KEY[i];
	model.setObjective( nv, GRB_MAXIMIZE );

	if ( time > 0 )
		model.set(GRB_DoubleParam_TimeLimit, time );

	model.optimize();

    if ( model.get( GRB_IntAttr_Status ) == GRB_TIME_LIMIT )
    {
		logger( __func__  + string( " : " ) + to_string( rounds ) + 
				string( " | EXPAND " ) );
        return EXPAND;
    }

	else
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
		if ( solCount >= MAX )
		{
			cerr << "solCount value  is too big !" << endl;
			exit(-1);
		}
		double xtime = model.get( GRB_DoubleAttr_Runtime );

		logger(__func__  + string( " : " ) + to_string( rounds ) + string( " | "
					) + to_string( xtime ) + string( " | "
					) + to_string( solCount ) );

		if ( solCount > 0 ) 
		{
			bitset<128> start;
			for ( int i = 0; i < solCount; i++ )
			{
				model.set(GRB_IntParam_SolutionNumber, i );

				for ( int j = 0; j < 128; j++ ) 
					if ( round( KEY[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
						start[j] = 1;
					else 
						start[j] = 0;
				counterMap[start]++;
			}
			return SOLUTION;
        }
		else
			return NOSOLUTION;

    }
}

/*
int main()
{
	vector<int> CUBE{ 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
		18, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
		37, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
		57, 58, 59, 60, 62, 63, 64, 65, 67, 68, 69, 70, 71, 74, 75, 76, 77, 79,
		80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
		98, 99, 100, 102, 103, 104, 105, 107, 108, 111, 112, 113, 114, 115, 116,
		117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127 };

	bitset<128> cube;
	for ( int i = 0; i < 128; i++ )
		if ( count( CUBE.begin(), CUBE.end(), i ) )
			cube[i] = 1;
	    else
			cube[i] = 0;

	map<bitset<128>, int, CMPS<128>> counterMap;

	kreyvium( 892, cube, counterMap );

	for( auto it : counterMap ) 
	{
		cout << it.first << " " << it.second << endl;
		if ( it.second % 2 == 1 )
		{
			for ( int i = 0; i < 128; i++ )
			{
				if ( it.first[i] == 1 )
					cout << "k" << i << " ";
			}
			cout << endl;
		}
	}

	vector<bitset<256>> terms;

	BackExpand( 100, terms );

	for ( auto it : terms )
	{
		for ( int i = 0; i < 256; i++ )
			if ( it[i] == 1 )
			{
				if ( i < 128 )
					cout << "k" << i << ' ';
				else
					cout << "v" << i - 128 << ' ';
			}
		cout << endl;
	}
}
*/










