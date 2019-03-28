//
// Created by Chris Broughall on 2019-03-28.
//


#include <cmath>
#include <memory>
#include "EulerMethod.hpp"
#include "RNG.hpp"


vector<double> Euler::genStockPrices (double S0, double T, double r, double sigma, int simSize)
{

    //double dt = 0.5/simSize;
    double dt = 1.0/365.0;
    double t;



    double St = S0;
    vector<double> All(simSize,0.0);
    double StOld = St;

    int j = 0;

        for (int i = 0; i < simSize; i++) {

            t = 0.0;
            St = S0;
            shared_ptr<RNG> rng=make_shared<RNG>();
            vector<double> e= rng->rngUsingBM(simSize*round(T/dt) + 1);
            j = 0;
            StOld = S0;

            do{

                St = StOld + (StOld * dt*r) + (sigma * sqrt(dt) * StOld * e[j]);
                StOld = St;


                t = t + dt;
                j++;

            }while(t <= T);

            //std::cout << "The St is: " << St <<"\n";

            All[i] = St;

        }



    return All;

}
