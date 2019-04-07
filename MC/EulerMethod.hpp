//
// Created by Chris Broughall on 2019-03-28.
//

#ifndef MC_EULERMETHOD_HPP
#define MC_EULERMETHOD_HPP

#include <vector>

using namespace std;

class Euler
{
public:
    vector<double> genStockPrices(double S0, double T, double r, double sigma, int simSize);
    double getStockPrice(const double S0, const double T, const double r, const double sigma, const double dt);
    void writeToFile (const string fNamePrefix, unique_ptr<double[]> v, const int size, const int k);
};



#endif //MC_EULERMETHOD_HPP
