#include<iostream>
#include<fstream>
#include<vector>
#include<bitset>
#include<algorithm>
#include<thread>
#include<mutex>
#include<future>
#include"deg.h"
#include"trivium.h"
#include"log.h"
#include"thread_pool.h"

using namespace std;
using namespace thread_pool;
const int THREAD_NUM = 64;

mutex stream_mutex;
mutex map_mutex;
mutex term_mutex;

void filterTerms( const bitset<80> & cube, int rounds, vector<bitset<288>> & terms )
{ 
	int size0 = terms.size();
    vector<bitset<288>> tmp( terms.begin(), terms.end() );
    terms.clear();
    for ( auto it : tmp )
    {
        auto d = computeDegree( cube, rounds, it );
        if ( d >= cube.count() )
			terms.push_back( it );
    }
	int size1 = terms.size();

	logger( __func__ + string(": ") + to_string( size0 ) + string( "\t" ) +
			to_string( size1 ) );
}

void filterMap( map<bitset<288>, int, CMPS<288>> & mp )
{ 
	int size0 = mp.size();
	map< bitset<288>, int, CMPS<288>> tmp( mp );
    mp.clear();

    for ( auto it : tmp )
    {
		if ( it.second % 2 == 1 )
			mp[it.first] = it.second;
    }
	int size1 = mp.size();

	logger( __func__ + string(": ") + to_string( size0 ) + string( "\t" ) +
			to_string( size1 ) );
}

void printSol( int rounds, const bitset<288> & vector, const map<bitset<80>,
		int, CMPS<80>> & solutions )  
{
	logger( __func__ + string(": " ) + to_string( rounds ) + string("\t") +
			vector.to_string()  );

	string path = string ( "TERM/" ) + to_string( rounds ) + string ( ".txt" ); 
	ofstream os;
	os.open( path, ios::out | ios::app );
	os << rounds << "\t" << vector << endl;
	for ( auto it : solutions )
		os << it.first << "\t" << it.second << "\n";
	os << endl;
	os.close();
}

void SolutionSearcherWorker(const bitset<288> & vec,  int ROUND, const bitset<80> & cube,
		vector<bitset<288>> & layerTerms, float time, int singlethread )
{
	map<bitset<80>, int, CMPS<80>> solutions;
	 
	auto status = MidSolutionCounter( ROUND, cube, vec,
				solutions, time, singlethread );

	if ( status == EXPAND )
    {
		lock_guard<mutex> guard( map_mutex );
        layerTerms.push_back( vec );
    }
	else if ( status == SOLUTION )
	{
		lock_guard<mutex> guard( stream_mutex );
		printSol( ROUND, vec, solutions );
	}
}

void ExpandWorker( const vector<bitset<288>> & layerTerm, const bitset<80> & cube, 
		int start, int end, int step, 
		map<bitset<288>, int, CMPS<288>> & layerMap, int singlethread )
{
	map<bitset<288>, int, CMPS<288>>  threadLayerMap;
	 
	int END = end <= layerTerm.size() ? end : layerTerm.size();

	vector<bitset<288>> terms;


	for ( int i = start; i < END; i++ )
	{
	    terms.clear();
		SecondBackExpandPolynomial( step, layerTerm[i], terms, singlethread  );

		for ( auto it : terms )
			threadLayerMap[it] += 1;
	}

	lock_guard<mutex> guard( term_mutex );
	for ( auto & it : threadLayerMap )
		layerMap[it.first] += it.second;
}


int main()
{
    // change CubeIndex for other cube
    vector<int> cubeIndex { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
            27, 28, 29, 30, 31, 32, 34, 36, 38, 40, 42, 45, 47, 49, 51, 53, 55, 57, 60, 62,
            64, 66, 68, 70, 72, 75, 77, 79 };

	int ROUND = 845;

	bitset<80> cube;
	for ( int i = 0; i < 80; i++ )
		if ( count( cubeIndex.begin(), cubeIndex.end(), i ) )
			cube[i] = 1;
	    else
			cube[i] = 0;

	int step = 200;

	float time = 0;
	map<bitset<80>, int, CMPS<80>> solutions;

    vector<bitset<288>> terms;
	vector<bitset<288>> layerTerms;
	map<bitset<288>, int, CMPS<288>> layerMap;

    BackExpandPolynomial( step, terms );

    ROUND = ROUND - step;

    filterTerms( cube, ROUND, terms );
	logger( __func__ + string( "First expand" ) );

	layerMap.clear();
	for ( auto it : terms )
	    layerMap[it] ++;

    ThreadPool thread_pool{};
	vector<bitset<288>> mapTmp;
	while ( true )
	{
	    // push all the terms into LayerMap
		layerTerms.clear();
		filterMap( layerMap );

		int singlethread = 2;

		mapTmp.clear();
		for ( auto& it : layerMap )
			mapTmp.push_back( it.first );

		filterTerms( cube, ROUND, mapTmp );

        if ( ROUND > 400 )
            time = 10;
        else if ( ROUND > 100 )
            time = 30;
        else if ( ROUND > 400 )
            time = 180;
        else if ( ROUND > 300 )
            time = 360;
        else if ( ROUND > 200 )
            time = 720;
        else if ( ROUND > 100 )
            time = 1200;
        else if ( ROUND > 20 )
            time = 3600;
        else
            time = 0;

        vector<std::future<void>> futures;
        for (auto & it : mapTmp ) 
            futures.emplace_back(thread_pool.Submit( SolutionSearcherWorker, it,
                        ROUND, ref(cube), ref( layerTerms ), time, singlethread ) );
        for (auto & it : futures)
            it.get();

		logger( string( "layerTermsSize " ) + to_string( layerTerms.size() ) );

		if ( layerTerms.size() == 0 )
			break;

		if ( ROUND <= 0 )
			cerr << "Error " << ROUND << endl;
		else if (ROUND < 10)
			time = 0;

        step = 0;
		layerMap.clear();

        vector<thread> threads;
        while( layerMap.size() < 100000 )
        {
            step ++;

		    layerMap.clear();

		    int size = layerTerms.size();

		    int base = ( size / THREAD_NUM ) >= 1 ? (size / THREAD_NUM) :
			1;
		    int size_multiple = base * THREAD_NUM;

		    threads.clear();
		    singlethread = 2;

		    for ( int i = 0; i < THREAD_NUM; i++ )
			    threads.push_back( thread ( ExpandWorker, ref(
							layerTerms ), ref ( cube ), base * i, base * (i + 1), step , 
							ref (layerMap ), singlethread ) );
		    for ( auto & th : threads ) th.join();

		    if ( size > size_multiple )
		    {
                threads.clear();
                for ( int i = 0; i < THREAD_NUM; i++ )
                    threads.push_back( thread ( ExpandWorker, ref(
                                layerTerms ), ref( cube ),  size_multiple + i, size_multiple + i + 1, step, 
                                ref( layerMap ), singlethread ) );
                for ( auto & th : threads ) th.join();
            }
        }

		ROUND = ROUND - step;
		logger( "Step New " + to_string( step ) );
		logger( "Round New " + to_string( ROUND ) );
		logger( string( "layerMapSize " ) + to_string( layerMap.size() ) );
	}

	cout << "Success!" << endl;
}
