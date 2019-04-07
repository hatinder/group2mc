//
// Created by hsingh9 on 23/03/2019.
//

#ifndef MC_MONTECARLO_HPP
#define MC_MONTECARLO_HPP


#include <vector>
using namespace std;

class MonteCarlo
{
public:
    double genStockPrice (double S0, double T, double r, double sigma);
    double genStockPrice2 (double S0, double T, double r, double sigma, std::vector<double> epsilon);
    vector<double> genStockPrices (double S0, double T, double r, double sigma,int simSize);
    double genStockPrices (double S0, double T, double r, double sigma, double dt, int simsize);
    void writeToFile (const string fNamePrefix, unique_ptr<double[]> v, const int size, const int k);
};


#endif //MC_MONTECARLO_HPP
