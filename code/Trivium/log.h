#ifndef __LOG_H__
#define __LOG_H__
#include<chrono>
#include<iostream>
#include<fstream>
#include<string>

using namespace std;

typedef bool TIMETYPE;

const int DATE = false;
const int TIME = true;

inline string getCurrentSystemTime(TIMETYPE type)
{
    auto tt = chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm* ptm = localtime(&tt);
    char date[60] = { 0 };
    if ( type == DATE )
        sprintf( date, "%d-%02d-%02d", (int)ptm->tm_year + 1900, (int)ptm->tm_mon + 1, (int)ptm->tm_mday );
    else 
        sprintf( date, "%02d:%02d:%02d", (int) ptm->tm_hour, (int)ptm->tm_min, (int)ptm -> tm_sec);

    return string(date);
}

inline void logger( string logMsg )
{
    string filePath = "LOG/log_" + getCurrentSystemTime( DATE ) + ".log"; 
    string now = getCurrentSystemTime(TIME);
    ofstream os ( filePath.c_str(), ios_base::out | ios_base::app ); 
    os << now << "\t" << logMsg << "\n";
    os.close();
}

#endif
