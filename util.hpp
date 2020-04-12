#ifndef _CPP_UTIL_
#define _CPP_UTIL_

#include<iostream>
#include<string>
#include<sstream>
//#include<boost/numeric/ublas/vector.hpp>

using namespace std;

//////////////////// utils/////////////////////////////////
template < class T> 
void convert_from_string(T& value, const string& s)
{
    stringstream ss(s);
    ss >> value;
}

double norm(const std::vector<double>& v1, const std::vector<double>& v2) {
    assert (v1.size() == v2.size());
    double sum = 0;
    for (size_t i=0; i<v1.size(); ++i)
    {
        double minus = v1[i] - v2[i];
        double r = minus * minus;
        sum += r;
    }

    return sqrt(sum);
}


double norm_1(const std::vector<double>& v1, const std::vector<double>& v2) {
    assert (v1.size() == v2.size());
    double sum = 0;
    for (size_t i=0; i<v1.size(); ++i)
    {
        double minus = abs(v1[i] - v2[i]);
        sum += minus;
    }

    return sum;
}

#endif
