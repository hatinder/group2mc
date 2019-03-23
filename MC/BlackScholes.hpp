//
// Created by hsingh9 on 23/03/2019.
//

#ifndef MC_BLACKSCHOLES_HPP
#define MC_BLACKSCHOLES_HPP

#include <vector>
using namespace std;
class BlackScholes
{
public:
    vector<double> calcExactValues (double S0, double K, double r, double sigma, double T);
    double calcExactDelta();
};


#endif //MC_BLACKSCHOLES_HPP
