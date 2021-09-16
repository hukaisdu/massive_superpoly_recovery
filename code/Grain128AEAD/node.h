#ifndef __NODE_H__
#define __NODE_H__
#include<vector>
#include<bitset>
#include<map>

using namespace std;

enum STATUS { SOLUTION, NOSOLUTION, EXPAND };

template<size_t N>
bool CMP( const bitset<N> & a, const bitset<N> & b )
{
        for ( int i = 0; i < N; i++ )
            if ( a[i] < b[i] ) return true;
            else if ( a[i] > b[i] ) return false;
        return false; // equal
}

template<size_t N>
struct CMPS
{
    bool operator()( const bitset<N> & a, const bitset<N> & b ) const
    {
        for ( int i = 0; i < N; i++ )
            if ( a[i] < b[i] ) return true;
            else if ( a[i] > b[i] ) return false;
        return false; // equal
    }
};

struct STAIRS
{
    int nextStep;
	float time;
};

struct Node
{
    int _rnd;
    STAIRS stair;    
    bitset<288> _vector;
    STATUS _status;
    map< bitset<80>, int, CMPS<80> > _solutions;
    vector<Node> _child;
};

struct cmpNode
{
    bool operator() ( const Node & a, const Node & b ) const 
    {
        if ( a._rnd < b._rnd )
            return true;
        else if ( a._rnd > b._rnd )
            return false;
        return CMP( a._vector, b._vector );
    }
};

#endif
