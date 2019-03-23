//
// Created by hsingh9 on 23/03/2019.
//

#include <iostream>
#include "BlackScholes.hpp"
#include "stats.hpp"

using namespace stats;
using namespace std;

vector<double> BlackScholes::calcExactValues (double S0, double K, double r, double sigma, double T)
{
    vector<double> values(2,0.0);
    double d1=(log(S0/K) + (r+(sigma*sigma)/2) * T)/(sigma*sqrt(T));
    double d2=(log(S0/K) + (r-(sigma*sigma)/2) * T)/(sigma*sqrt(T));
    double Nd1=pnorm(d1,0.0,1.0);
    double Nd2=pnorm(d2,0.0,1.0);
    values[0]=S0*Nd1-K*exp(-r*T)*Nd2;
//    cout<<"d1: "<<d1<<", d2:"<<d2<<", Nd1:"<<Nd1<<", Nd2:"<<Nd2<<endl;
    values[1]=Nd1;
    return values;
}

