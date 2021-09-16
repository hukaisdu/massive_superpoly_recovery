#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"

const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

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

int SecondBackExpandPolynomial( int rounds, const bitset<288> & last, vector<bitset<288> > & term, int threads  )
{
	if ( term.size() > 0 ) 
	{
		cerr <<  __func__ << " Term is not empty" << endl;
		exit(-1);
	}

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2); 
    env.set(GRB_IntParam_Threads, threads );
    env.set(GRB_IntParam_PoolSolutions, MAX); 

    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> works = s;
    for (int r = 0; r < rounds; r++) 
    {
        triviumCore(model, works, 65, 170, 90, 91, 92);
        triviumCore(model, works, 161, 263, 174, 175, 176);
        triviumCore(model, works, 242, 68, 285, 286, 287);
            
        vector<GRBVar> temp = works;
        for (int i = 0; i < 288; i++) 
            works[(i + 1) % 288] = temp[i];
    }

    // output constraints
    for ( int i = 0; i < 288; i++ )
        if ( last[i] == 1 )
            model.addConstr( works[i] == 1 );
        else
            model.addConstr( works[i] == 0 );

	//logger(__func__ );
    model.optimize();

    map<bitset<288>, int, CMPS<288>> counterMap; 

    if( model.get( GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
		if ( solCount >= MAX )
		{
			cerr << "solCount value  is too big !" << endl;
			exit(-1);
		}
        bitset<288> start;
        for ( int i = 0; i < solCount; i++ )
        {
            model.set(GRB_IntParam_SolutionNumber, i );

            for ( int j = 0; j < 288; j++ ) 
                if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;
            counterMap[start]++;
        }
    }
    else if( model.get( GRB_IntAttr_Status ) == GRB_INFEASIBLE )
    {
        exit(-2);
    }
    else
    {
        exit(-1);
    }

    for ( auto it : counterMap )
        if ( it.second % 2 == 1 )
            term.push_back( it.first );
    logger( string( __func__ ) + string( ":" )  + string( "SecondExpand Terms: " ) + to_string( term.size() ) );
    return 0;
}

int BackExpandPolynomial( int rounds, vector<bitset<288> > & term )
{
	if ( term.size() > 0 ) 
	{
		cerr << __func__ << " Term is not empty" << endl;
		exit(-1);
	}
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2); 
    env.set(GRB_IntParam_PoolSolutions, MAX); 

    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> works = s;
    for (int r = 0; r < rounds; r++) 
    {
        triviumCore(model, works, 65, 170, 90, 91, 92);
        triviumCore(model, works, 161, 263, 174, 175, 176);
        triviumCore(model, works, 242, 68, 285, 286, 287);
            
        vector<GRBVar> temp = works;
        for (int i = 0; i < 288; i++) 
            works[(i + 1) % 288] = temp[i];
    }

    // output constraints
    GRBLinExpr nk = 0;
    for ( int i = 0; i < 288; i++ )
        if ( (i == 65) || (i == 92) || (i == 161) || (i == 176) || (i == 242) || (i == 287))
            nk += works[i];
        else 
            model.addConstr( works[i] == 0);
    model.addConstr( nk == 1 );

    model.update();
	logger(__func__ );
    model.optimize();

    map<bitset<288>, int, CMPS<288>> counterMap; 

    if( model.get( GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
		if ( solCount >= MAX )
		{
			cerr << "solCount value  is too big !" << endl;
			exit(-1);
		}
        bitset<288> start;
        for ( int i = 0; i < solCount; i++ )
        {
            model.set(GRB_IntParam_SolutionNumber, i );

            for ( int j = 0; j < 288; j++ ) 
                if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;
            counterMap[start]++;
        }
    }
    else if( model.get( GRB_IntAttr_Status ) == GRB_INFEASIBLE )
    {
        exit(-2);
    }
    else
    {
        exit(-1);
    }

    for ( auto it : counterMap )
        if ( it.second % 2 == 1 )
            term.push_back( it.first );
    logger( string( __func__ ) + string( ":" )  + string( "Expand Terms: " ) + to_string( term.size() ) );
    return 0;
}

STATUS MidSolutionCounter( int rounds, const bitset<80> & cube, const
		bitset<288> & last, map<bitset<80>, int, CMPS<80>> & counterMap, float
		time, int thread )
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread );
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_MIPFocus, 3 );
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for ( int i = 0; i < 80; i++ )
        if ( cube[i] == 0  )
            model.addConstr( s[i + 93] == 0 );
        else
            model.addConstr( s[i + 93] == 1 );

    // key, last three bits
    for ( int i = 80; i < 93; i++ )
        model.addConstr( s[i] == 0 );
    for ( int i = 93 + 80; i < 285; i++ )
        model.addConstr( s[i] == 0 );

    vector<GRBVar> works = s;
    for (int r = 0; r < rounds; r++) 
    {
        triviumCore(model, works, 65, 170, 90, 91, 92);
        triviumCore(model, works, 161, 263, 174, 175, 176);
        triviumCore(model, works, 242, 68, 285, 286, 287);
            
        vector<GRBVar> temp = works;
        for (int i = 0; i < 288; i++) 
            works[(i + 1) % 288] = temp[i];
    }

    for ( int i = 0; i < 288; i++ )
        if ( last[i] == 1)
            model.addConstr( works[i] == 1 );
        else
            model.addConstr( works[i] == 0 );

    GRBLinExpr nk = 0;
    for ( int i = 0; i < 80; i++ )
        nk += s[i];
    model.setObjective( nk, GRB_MAXIMIZE );

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

		double time = model.get( GRB_DoubleAttr_Runtime );
		logger(__func__  + string( " : " ) + to_string( rounds ) + string( " | "
					) + to_string( time ) + string( " | "
					) + to_string( solCount ) );

		
        if ( solCount > 0 )
        {
            bitset<80> start;
            for ( int i = 0; i < solCount; i++ )
            {
                model.set(GRB_IntParam_SolutionNumber, i );
                for ( int j = 0; j < 80; j++ ) 
                    if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
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
