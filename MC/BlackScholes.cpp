//
// Created by hsingh9 on 23/03/2019.
//

#include <iostream>
#include "BlackScholes.hpp"
// Fix stats.hpp , download using command git clone -b master --single-branch https://github.com/kthohr/stats ./stats
// Download to folder : ../ i.e. one folder above using above given command.
#include "stats.hpp"
// Download this file to following folder :

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
    cout<<"Black Scholes European Call Price: "<<values[0]<<", Delta: "<<values[1]<<endl;
    return values;
}

vector<double> BlackScholes::calcAssetOrNothingExactValues (double S0, double K, double r, double sigma, double T)
{
    vector<double> values(2,0.0);
    double d1=(log(S0/K) + (r+(sigma*sigma)/2) * T)/(sigma*sqrt(T));
    double Nd1=pnorm(d1,0.0,1.0);
//    values[0]=S0*exp((-r)*T)*Nd1;
    values[0]=S0*Nd1;
//    values[1]=exp((-r)*T)*(Nd1 + exp(((-1)*S0*S0)/2)/(sigma*sqrt(2*M_PI*T)));
//    values[1]=Nd1 + exp(((-1)*d1*d1)/2)/(sigma*sqrt(2*M_PI*T));
    values[1]=Nd1 + dnorm(d1,0.0,1.0)/(sigma*sqrt(T));
    cout<<"Black Scholes Asset-Or-Nothing Call Price: "<<values[0]<<", Delta: "<<values[1]<<endl;
    return values;
}