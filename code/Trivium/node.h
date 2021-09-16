#ifndef __NODE_H__
#define __NODE_H__
#include<vector>
#include<bitset>
#include<map>

using namespace std;

enum STATUS { SOLUTION, NOSOLUTION, EXPAND };

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

#endif
