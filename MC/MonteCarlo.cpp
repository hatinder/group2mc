//
// Created by hsingh9 on 23/03/2019.
//

#include <cmath>
#include <memory>
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
        AllST[j]=St;
    }
    return AllST;
}
