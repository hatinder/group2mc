//
// Created by hsingh9 on 23/03/2019.
//

#include <cmath>
#include <memory>
#include <sstream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include "MonteCarlo.hpp"
#include "RNG.hpp"

double MonteCarlo::genStockPrice (double S0, double T, double r, double sigma)
{
    double dt=1.0/365.0,t;
    double StNext,St;
    St=S0;
    t=0.0;
    do{
        shared_ptr<RNG> rng=make_shared<RNG>();
        vector<double> e= rng->rngUsingMTE(0);
        int i=0;
        do{
            double epsilon=e[i];
            St= St+St*exp( ( (r-(1.0/2.0) * sigma * sigma) )*dt + sigma*sqrt(dt)*epsilon );
            t+=dt;
            i++;
        }while(t<T && i<e.size());
    }while(t<T);
    return StNext;
}

double MonteCarlo::genStockPrice2 (double S0, double T, double r, double sigma, vector<double> epsilon)
{
    double dt=1.0/365.0,t,dtsqrt;
    dtsqrt=sqrt(dt);
    double StNext,St;
    St=S0;
    t=0.0;
    do{
        int i=0;
        do{
            double a=(r-(1.0/2.0) * sigma * sigma);
            double b=sigma*dtsqrt*epsilon[i];
            St=St*exp( a*dt + b );
            cout<<"a: "<<a<<", b: "<<b<<", epsilon: "<<epsilon[i]<<", St: "<<St<<endl;
            t+=dt;
            i++;
        }while(t<T && i<epsilon.size());
        cout<<"i: "<<i<<endl;
    }while(t<T);

    return St;
}

vector<double> MonteCarlo::genStockPrices (double S0, double T, double r, double sigma,int simSize)
{
    double dt=1.0/365.0,t;
    double St;
    double a=(r-(1.0/2.0) * sigma * sigma)*dt;
    double b=sigma*sqrt(dt);
    vector<double> AllST(simSize,0.0);
    shared_ptr<RNG> rng=make_shared<RNG>();
//    vector<double> e= rng->rngUsingMTE(simSize*round(T/dt) + 1);
    vector<double> e= rng->rngUsingBM(simSize*round(T/dt) + 1);
    int i=0;
    for (int j = 0; j < simSize; ++j)
    {
        St=S0;
        t=0.0;
        do{
            do{
                St= St*exp( a + b*e[i] );
                t+=dt;
                i++;
            }while(t<T && i<e.size());
        }while(t<T);
//        cout<<St<<endl;
        AllST[j]=St;
    }
    return AllST;
}

double MonteCarlo::genStockPrices (double S0, double T, double r, double sigma, double dt, int simsize)
{
    double St=S0,t=0.0;
    double a=(r-(1.0/2.0) * sigma * sigma)*dt;
    double b=sigma*sqrt(dt);
    unique_ptr<double[]> bmVal{new double[2]()};
    unique_ptr<RNG> rng;
    while(t<=T)
    {
        bmVal= rng->rngUsingStatsBM(0);
        St=St*exp(a+b*bmVal[0]);
        St=St*exp(a+b*bmVal[1]);
        t=t+2*dt;
    }
//    cout<<"size="<<size<<endl;
    return St;
}

void MonteCarlo::writeToFile (const string fNamePrefix, unique_ptr<double[]> v, const int size, const int k)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    oFileStream<<setw(12)<<"x"<<setw(12)<<"y"<<endl;
    oFileStream<<endl;
    for (int i = 0; i < size; ++i)
    {
        oFileStream<<setw(12)<<i<<setw(12)<<v[i]<<endl;
    }
    oFileStream.close();

}

