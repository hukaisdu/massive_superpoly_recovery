#include<iostream>
#include<fstream>
#include<vector>
#include<bitset>
#include<algorithm>
#include<thread>
#include<mutex>
#include<future>
#include"log.h"
#include"thread_pool.h"
#include"grain.h"

using namespace std;
using namespace thread_pool;
const int THREAD_NUM = 64;
mutex stream_mutex;
mutex map_mutex;
mutex term_mutex;

void filterMap( map<bitset<256>, int, CMPS<256>> & mp )
{ 
	int size0 = mp.size();
	map< bitset<256>, int, CMPS<256>> tmp( mp );
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

void printSol( int rounds, const bitset<256> & vector, const map<bitset<128>,
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


void SolutionSearcherWorker(const bitset<256> & vec,  int ROUND, const bitset<96> & cube,
		vector<bitset<256>> & layerTerms, float time, int singlethread )
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

void ExpandWorker( const vector<bitset<256>> & layerTerm, const bitset<96> & cube, 
		int start, int end, int step, 
		map<bitset<256>, int, CMPS<256>> & layerMap, int singlethread )
{
	//map<bitset<80>, int, CMPS<80>> solutions;
	map<bitset<256>, int, CMPS<256>>  threadLayerMap;
	 
	int END = end <= layerTerm.size() ? end : layerTerm.size();

	vector<bitset<256>> terms;


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


int main()
{
	bitset<96> cube;
	for ( int i = 0; i < 96; i++ )
		cube[i] = 1;

	int ROUND = 192;
	int step = 31;

	float time = 0;
	map<bitset<128>, int, CMPS<128>> solutions;

    vector<bitset<256>> terms;

	vector<bitset<256>> layerTerms;

	map<bitset<256>, int, CMPS<256>> layerMap;

    BackExpandPolynomial( step, terms );

    ROUND = ROUND - step;

    //filterTerms( cube, ROUND, terms );
	logger( __func__ + string( "First expand" ) );

	layerMap.clear();

	for ( auto it : terms )
	    layerMap[it] ++;

    ThreadPool thread_pool{};

	vector<bitset<256>> mapTmp;

	while ( true )
	{
	// push all the terms into LayerMap
		layerTerms.clear();
		filterMap( layerMap );

		int singlethread = 2;

		mapTmp.clear();

		for ( auto& it : layerMap )
			mapTmp.push_back( it.first );

		ofstream f; 
		f.open( string ( "STORE/store_" ) + to_string( ROUND ), ios::out );
		f << ROUND << endl;
		f << mapTmp.size() << endl;
		for ( auto it : mapTmp )
			f << it << endl;
		f.close();

        if ( ROUND > 120 )
            time = 60;
        else if ( ROUND > 110 )
            time = 120;
        else if ( ROUND > 100 )
            time = 180;
        else if ( ROUND > 33 )
            time = 360;
        else if ( ROUND > 10 )
            time = 360;
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
		//}
		logger( string( "layerMapSize " ) + to_string( layerMap.size() ) );
	}

	cout << "Success!" << endl;
}
