//
// Created by Chris Broughall on 2019-03-28.
//


#include <cmath>
#include <memory>
#include "EulerMethod.hpp"
#include "RNG.hpp"


vector<double> Euler::genStockPrices (double S0, double T, double r, double sigma, int simSize)
{

    double dt = 0.5/simSize;


    shared_ptr<RNG> rng=make_shared<RNG>();
    vector<double> e= rng->rngUsingBM(simSize);
    double St = S0;
    vector<double> All(simSize,0.0);
    for (int i = 0; i < simSize; i++)
    {
        St = St + St*dt + sigma*sqrt(dt)*e[i];
        All[i] = St;
    }

    return All;

}
