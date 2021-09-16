#include<iostream>
#include<fstream>
#include<vector>
#include<bitset>
#include<algorithm>
#include<thread>
#include<mutex>
#include<future>
#include"log.h"
#include"deg.h"
#include"thread_pool.h"
#include"kreyvium.h"

using namespace std;
using namespace thread_pool;
mutex stream_mutex;
mutex map_mutex;
mutex term_mutex;
int THREAD_NUM = 64;

void filterMap( map<bitset<544>, int, CMPS<544>> & mp )
{ 
	int size0 = mp.size();
	map< bitset<544>, int, CMPS<544>> tmp( mp );
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

void printSol( int rounds, const bitset<544> & vector, const map<bitset<128>,
		int, CMPS<128>> & solutions )  
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

void SolutionSearcherWorker( const bitset<544> & vec,  int ROUND, 
		const bitset<128> & cube,
		vector< bitset<544> > & layerTerms, float time, int singlethread )
{
	map<bitset<128>, int, CMPS<128>> solutions;
	 
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

void ExpandWorker( const vector<bitset<544>> & layerTerm, const bitset<128> & cube, 
		int start, int end, int step, 
		map<bitset<544>, int, CMPS<544>> & layerMap, int singlethread )
{
	//map<bitset<80>, int, CMPS<80>> solutions;
	map<bitset<544>, int, CMPS<544>>  threadLayerMap;
	 
	int END = end <= layerTerm.size() ? end : layerTerm.size();

	vector<bitset<544>> terms;


	for ( int i = start; i < END; i++ )
	{
		//MidSolutionCounter( ROUND, cube, layerMap[i],
	//			solutions, time, singlethread ); 
	    terms.clear();
		SecondBackExpandPolynomial( step, layerTerm[i], terms, singlethread  );

		//filterTerms( cube, ROUND, terms );

		for ( auto it : terms )
			threadLayerMap[it] += 1;
	}

	lock_guard<mutex> guard( term_mutex );
	for ( auto & it : threadLayerMap )
		layerMap[it.first] += it.second;
}

void filterTerms( const bitset<128> & cube, int rounds, vector<bitset<544>> & terms )
{ 
	int size0 = terms.size();
    vector<bitset<544>> tmp( terms.begin(), terms.end() );
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


int main()
{
	/* change cubeIndex to test other cube */
	vector<int> cubeIndex{ 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
		18, 19, 20, 21, 22,/**/  23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
		37, 38, /**/ 39, 40, 41, 42, 43, /**/ 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
		57, 58, 59, 60, 61, /**/ 62, 63, 64, 65, 67, 68, 69, 70, 71, 74, 75, 76, 77, 79,
		80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
		98, 99, 100, 102, 103, 104, 105, 107, 108, 111, 112, 113, 114, 115, 116,
		117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127 };

	bitset<128> cube{ 0 };
	for ( auto it : cubeIndex )
		cube[it] = 1;

	/*change the ROUND to test the 894 round */
	int ROUND = 893;
	int step = 200;

	float time = 0;
	map<bitset<128>, int, CMPS<128>> solutions;

    vector<bitset<544>> terms;

	vector<bitset<544>> layerTerms;

	map<bitset<544>, int, CMPS<544>> layerMap;

    BackExpandPolynomial( step, terms );

    ROUND = ROUND - step;

	/* use numeric mapping to filter out some terms that have no contribution to
	 * the superpoly
	 */
    filterTerms( cube, ROUND, terms );

	logger( __func__ + string( "First expand" ) );

	layerMap.clear();

	for ( auto it : terms )
	    layerMap[it] ++;

    ThreadPool thread_pool{};

	vector<bitset<544>> mapTmp;

	while ( true )
	{
		// push all the terms into LayerMap
		layerTerms.clear();
		filterMap( layerMap );

		int singlethread = 2;
		mapTmp.clear();
		for ( auto& it : layerMap )
			mapTmp.push_back( it.first );

		/* numeric mapping to filter out some terms */
        filterTerms( cube, ROUND, mapTmp );

		/* ChooseTi, change it to adjust the costs of time and memory */
        if ( ROUND > 600 )
            time = 60;
        else if ( ROUND > 500 )
            time = 120;
        else if ( ROUND > 400 )
            time = 180;
        else if ( ROUND > 300 )
            time = 360;
        else if ( ROUND > 200 )
            time = 720;
        else
            time = 0;

		/* compute small superpoly */
        vector<std::future<void>> futures;
        for (auto & it : mapTmp ) 
            futures.emplace_back( thread_pool.Submit( SolutionSearcherWorker, it,
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
		/* ChooseRi, change 10000 to adjust the memory and time costs */
        while( layerMap.size() < 10000 )
        {
            step ++;
		    layerMap.clear();

		    int size = layerTerms.size();

		    int base = ( size / THREAD_NUM ) >= 1 ? (size / THREAD_NUM) :
			1;
		    int size_multiple = base * THREAD_NUM;

		    threads.clear();
		    singlethread = 2;

			// expand the terms 
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
