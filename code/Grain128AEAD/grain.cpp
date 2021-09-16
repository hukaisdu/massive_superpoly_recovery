#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"

using namespace std;
const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

GRBVar tap(GRBModel& model, GRBVar& x) {

  GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar tmp[2] = { y,z };
  model.addGenConstrOr(x, tmp, 2);
  x = z;
  return y;
}

GRBVar funcH(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s ) 
{
  GRBVar b12x = tap(model, b[12]);
  GRBVar s8 = tap(model, s[8]);
  model.addConstr(b12x == s8);

  GRBVar s13 = tap(model, s[13]);
  GRBVar s20 = tap(model, s[20]);
  model.addConstr(s13 == s20);

  GRBVar b95x = tap(model, b[95]);
  GRBVar s42 = tap(model, s[42]);
  model.addConstr(b95x == s42);

  GRBVar s60 = tap(model, s[60]);
  GRBVar s79 = tap(model, s[79]);
  model.addConstr(s60 == s79);

  GRBVar b12y = tap(model, b[12]);
  GRBVar b95y = tap(model, b[95]);
  GRBVar s94 = tap(model, s[94]);
  model.addConstr(b12y == b95y);
  model.addConstr(b12y == s94);

  GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(y == b12x + s13 + b95x + s60 + b12y);

  return y;
}

GRBVar funcO(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s ) {
  GRBVar s93 = tap(model, s[93]);
  GRBVar b2 = tap(model, b[2]);
  GRBVar b15 = tap(model, b[15]);
  GRBVar b36 = tap(model, b[36]);
  GRBVar b45 = tap(model, b[45]);
  GRBVar b64 = tap(model, b[64]);
  GRBVar b73 = tap(model, b[73]);
  GRBVar b89 = tap(model, b[89]);

  GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(y == s93 + b2 + b15 + b36 + b45 + b64 + b73 + b89);  

  return y;
}

GRBVar funcF(GRBModel& model, vector<GRBVar>& s) {
  GRBVar s0 = tap(model, s[0]);
  GRBVar s7 = tap(model, s[7]);
  GRBVar s38 = tap(model, s[38]);
  GRBVar s70 = tap(model, s[70]);
  GRBVar s81 = tap(model, s[81]);
  GRBVar s96 = tap(model, s[96]);

  GRBVar f = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(f == s0 + s7 + s38 + s70 + s81 + s96);

  return f;
}

GRBVar funcG(GRBModel& model, vector<GRBVar>& b) 
{
  // nonlinear
  GRBVar b26 = tap(model, b[26]);
  GRBVar b56 = tap(model, b[56]);
  GRBVar b91 = tap(model, b[91]);
  GRBVar b96 = tap(model, b[96]);

  GRBVar b3 = tap(model, b[3]);
  GRBVar b67 =  tap(model, b[67]);
  model.addConstr(b3 == b67);

  GRBVar b11 = tap(model, b[11]);
  GRBVar b13 = tap(model, b[13]);
  model.addConstr(b11 == b13);

  GRBVar b17 = tap(model, b[17]);
  GRBVar b18 = tap(model, b[18]);
  model.addConstr(b17 == b18);

  GRBVar b27=  tap(model, b[27]);
  GRBVar b59= tap(model, b[59]);
  model.addConstr(b27 == b59);

  GRBVar b40=  tap(model, b[40]);
  GRBVar b48=  tap(model, b[48]);
  model.addConstr(b40 == b48);

  GRBVar b61=  tap(model, b[61]);
  GRBVar b65=  tap(model, b[65]);
  model.addConstr(b61 == b65);

  GRBVar b68=  tap(model, b[68]);
  GRBVar b84=  tap(model, b[84]);
  model.addConstr(b68 == b84);

  GRBVar b88=  tap(model, b[88]);
  GRBVar b92=  tap(model, b[92]);
  GRBVar b93=  tap(model, b[93]);
  GRBVar b95=  tap(model, b[95]);
  model.addConstr(b88 == b92);
  model.addConstr(b88 == b93);
  model.addConstr(b88 == b95);

  GRBVar b22=  tap(model, b[22]);
  GRBVar b24=  tap(model, b[24]);
  GRBVar b25=  tap(model, b[25]);
  model.addConstr(b22 == b24);
  model.addConstr(b22 == b25);

  GRBVar b70=  tap(model, b[70]);
  GRBVar b78=  tap(model, b[78]);
  GRBVar b82=  tap(model, b[82]);
  model.addConstr(b70 == b78);
  model.addConstr(b70 == b82);

  // nonlinear feed back
  GRBVar g = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(g == b[0] + b26 + b56 + b91 + b96 + b3 + b11 + b17 + b27 + b40 + b61 + b68 + b88 + b22 + b70);

  return g;
}


int grain( int rounds, const bitset<96>& cube )
{
	//gurobi
    // Create the environment
    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, 2000000000 );
    //env.set(GRB_IntParam_PoolSolutions, 1 );
    env.set(GRB_IntParam_Threads, 8);

    // Create the model
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> s(128);
    vector<GRBVar> b(128);
    for (int i = 0; i < 128; i++)
	{
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
        b[i] = model.addVar(0, 1, 0, GRB_BINARY);
    }

    // IV constraint
    for (int i = 0; i < 96; i++)
        if (cube[i] == 1)
      	    model.addConstr(s[i] == 1);
        else
			model.addConstr(s[i] == 0);

    model.addConstr(s[127] == 0);

    // Round function

	vector<GRBVar> workb = b;
	vector<GRBVar> works = s;
	works = s;
    for (int r = 0; r <= rounds; r++) 
	{
        if (r < rounds) 
		{
			GRBVar s0 = works[0];
            GRBVar h = funcH(model, workb, works);
            GRBVar o = funcO(model, workb, works);

            GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
            model.addConstr(z == h + o);

            GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
            GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
            GRBVar tmpVar[2] = { z1, z2 };
            model.addGenConstrOr(z, tmpVar, 2);

            GRBVar f = funcF(model, works);
            GRBVar g = funcG(model, workb);

            GRBVar news = model.addVar(0, 1, 0, GRB_BINARY);
            model.addConstr(news == z1 + f);

            GRBVar newb = model.addVar(0, 1, 0, GRB_BINARY);
            model.addConstr(newb == z2 + g + works[0]);

			vector<GRBVar> tmpb = workb;
			vector<GRBVar> tmps = works;

            for (int i = 0; i < 127; i++) 
			{
                workb[i] = tmpb[i + 1];
                works[i] = tmps[i + 1];
            }
            workb[127] = newb;
            works[127] = news;

            // remove (s[r][0], z[r]) = (1,1) 
            model.addConstr((1 - s0 ) + (1 - z) >= 1);
        }
        else
		{
            GRBVar h = funcH(model, workb, works );
            GRBVar o = funcO(model, workb, works );

            GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
            model.addConstr(z == h + o);

            model.addConstr(z == 1);

            for (int i = 0; i < 128; i++) 
			{
                model.addConstr(workb[i] == 0);
                model.addConstr(works[i] == 0);
            }
        }
    }

    GRBLinExpr sumc = 0;
    for (int i = 0; i < 128; i++) {
        sumc += workb[i];
    }
    model.setObjective(sumc, GRB_MAXIMIZE);

    //
    model.update();

    model.optimize();

    if ( model.get(GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
        return solCount ;
    }
	else
	{
		cout << "Infeasible" << endl;
		return model.get( GRB_IntAttr_Status );
	}
}

STATUS MidSolutionCounter( int rounds, const bitset<96> & cube, const
		bitset<256> & last, map<bitset<128>, int, CMPS<128>> & counterMap, float
		time, int thread )
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread );
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_MIPFocus, 3 );
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> s(128);
    vector<GRBVar> b(128);
    for (int i = 0; i < 128; i++)
	{
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
        b[i] = model.addVar(0, 1, 0, GRB_BINARY);
    }

    // IV constraint
    for (int i = 0; i < 96; i++)
        if (cube[i] == 1)
      	    model.addConstr(s[i] == 1);
        else
			model.addConstr(s[i] == 0);

    model.addConstr(s[127] == 0);

    // Round function

	vector<GRBVar> workb = b;
	vector<GRBVar> works = s;
	works = s;
    for (int r = 0; r < rounds; r++) 
	{
		GRBVar s0 = works[0];
		GRBVar h = funcH(model, workb, works);
		GRBVar o = funcO(model, workb, works);

		GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(z == h + o);

		GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar tmpVar[2] = { z1, z2 };
		model.addGenConstrOr(z, tmpVar, 2);

		GRBVar f = funcF(model, works);
		GRBVar g = funcG(model, workb);

		GRBVar news = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(news == z1 + f);

		GRBVar newb = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(newb == z2 + g + works[0]);

		vector<GRBVar> tmpb = workb;
		vector<GRBVar> tmps = works;

		for (int i = 0; i < 127; i++) 
		{
			workb[i] = tmpb[i + 1];
			works[i] = tmps[i + 1];
		}
		workb[127] = newb;
		works[127] = news;

		// remove (s[r][0], z[r]) = (1,1) 
		model.addConstr((1 - s0 ) + (1 - z) >= 1);
    }
    // output constraints
    for ( int i = 0; i < 128; i++ )
        if ( last[i] == 1 )
            model.addConstr( workb[i] == 1 );
        else
            model.addConstr( workb[i] == 0 );

    for ( int i = 0; i < 128; i++ )
        if ( last[i + 128] == 1 )
            model.addConstr( works[i] == 1 );
        else
            model.addConstr( works[i] == 0 );

    GRBLinExpr sumc = 0;
    for (int i = 0; i < 128; i++) {
        sumc += b[i];
    }

    model.setObjective(sumc, GRB_MAXIMIZE);

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
            bitset<128> start;
            for ( int i = 0; i < solCount; i++ )
            {
                model.set(GRB_IntParam_SolutionNumber, i );
                for ( int j = 0; j < 128; j++ ) 
                    if ( round( b[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
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

int BackExpandPolynomial( int rounds,
		vector<bitset<256> > & term  )
{
	if ( term.size() > 0 ) 
	{
		cerr <<  __func__ << " Term is not empty" << endl;
		exit(-1);
	}

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2); 
    env.set(GRB_IntParam_PoolSolutions, MAX); 

    GRBModel model = GRBModel(env);


    // Create variables
    vector<GRBVar> s(128);
    vector<GRBVar> b(128);
    for (int i = 0; i < 128; i++)
	{
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
        b[i] = model.addVar(0, 1, 0, GRB_BINARY);
    }

    // Round function
	vector<GRBVar> workb = b;
	vector<GRBVar> works = s;
	works = s;
    for (int r = 0; r < rounds; r++) 
	{
		GRBVar s0 = works[0];
		GRBVar h = funcH(model, workb, works);
		GRBVar o = funcO(model, workb, works);

		GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(z == h + o);

		GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar tmpVar[2] = { z1, z2 };
		model.addGenConstrOr(z, tmpVar, 2);

		GRBVar f = funcF(model, works);
		GRBVar g = funcG(model, workb);

		GRBVar news = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(news == z1 + f);

		GRBVar newb = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(newb == z2 + g + works[0]);

		vector<GRBVar> tmpb = workb;
		vector<GRBVar> tmps = works;

		for (int i = 0; i < 127; i++) 
		{
			workb[i] = tmpb[i + 1];
			works[i] = tmps[i + 1];
		}
		workb[127] = newb;
		works[127] = news;

		// remove (s[r][0], z[r]) = (1,1) 
		model.addConstr((1 - s0 ) + (1 - z) >= 1);
    }

	// output constraints
	GRBVar h = funcH(model, workb, works );
	GRBVar o = funcO(model, workb, works );

	GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
	model.addConstr(z == h + o);

	model.addConstr(z == 1);

	for (int i = 0; i < 128; i++) 
	{
		model.addConstr(workb[i] == 0);
		model.addConstr(works[i] == 0);
	}


	//logger(__func__ );
    model.optimize();

    map<bitset<256>, int, CMPS<256>> counterMap; 

    if( model.get( GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
		if ( solCount >= MAX )
		{
			cerr << "solCount value  is too big !" << endl;
			exit(-1);
		}
        bitset<256> start;
        for ( int i = 0; i < solCount; i++ )
        {
            model.set(GRB_IntParam_SolutionNumber, i );

            for ( int j = 0; j < 128; j++ ) 
                if ( round( b[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;

            for ( int j = 0; j < 128; j++ ) 
                if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j + 128] = 1;
                else 
                    start[j + 128] = 0;
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

int SecondBackExpandPolynomial( int rounds, const bitset<256> & last,
		vector<bitset<256> > & term, int threads  )
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

    // Create variables
    vector<GRBVar> s(128);
    vector<GRBVar> b(128);
    for (int i = 0; i < 128; i++)
	{
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
        b[i] = model.addVar(0, 1, 0, GRB_BINARY);
    }

    // Round function

	vector<GRBVar> workb = b;
	vector<GRBVar> works = s;
	works = s;
    for (int r = 0; r < rounds; r++) 
	{
		GRBVar s0 = works[0];
		GRBVar h = funcH(model, workb, works);
		GRBVar o = funcO(model, workb, works);

		GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(z == h + o);

		GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar tmpVar[2] = { z1, z2 };
		model.addGenConstrOr(z, tmpVar, 2);

		GRBVar f = funcF(model, works);
		GRBVar g = funcG(model, workb);

		GRBVar news = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(news == z1 + f);

		GRBVar newb = model.addVar(0, 1, 0, GRB_BINARY);
		model.addConstr(newb == z2 + g + works[0]);

		vector<GRBVar> tmpb = workb;
		vector<GRBVar> tmps = works;

		for (int i = 0; i < 127; i++) 
		{
			workb[i] = tmpb[i + 1];
			works[i] = tmps[i + 1];
		}
		workb[127] = newb;
		works[127] = news;

		// remove (s[r][0], z[r]) = (1,1) 
		model.addConstr((1 - s0 ) + (1 - z) >= 1);
    }

    // output constraints
    for ( int i = 0; i < 128; i++ )
        if ( last[i] == 1 )
            model.addConstr( workb[i] == 1 );
        else
            model.addConstr( workb[i] == 0 );

    for ( int i = 0; i < 128; i++ )
        if ( last[i + 128] == 1 )
            model.addConstr( works[i] == 1 );
        else
            model.addConstr( works[i] == 0 );

	//logger(__func__ );
    model.optimize();

    map<bitset<256>, int, CMPS<256>> counterMap; 

    if( model.get( GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
		if ( solCount >= MAX )
		{
			cerr << "solCount value  is too big !" << endl;
			exit(-1);
		}
        bitset<256> start;
        for ( int i = 0; i < solCount; i++ )
        {
            model.set(GRB_IntParam_SolutionNumber, i );

            for ( int j = 0; j < 128; j++ ) 
                if ( round( b[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;

            for ( int j = 0; j < 128; j++ ) 
                if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j + 128] = 1;
                else 
                    start[j + 128] = 0;
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

/*
int main()
{
	vector<bitset<224>> term;
	BackExpandPolynomial( 10, term, 2 );
	for ( auto & it : term )
	{
		if ( it.count() == 0 )
			cout << "1" << endl;
		else
		{
			for ( int i = 0; i < 128; i++ )
				if ( it[i] == 1 )
					cout << "k" << i << ' ';
			for ( int i = 0; i < 96; i++ )
				if ( it[i + 128] == 1 )
					cout << "v" << i << ' ';
			cout << endl;
		}
	}
}
*/


