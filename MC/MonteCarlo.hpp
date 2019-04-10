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
    vector<double> genStockPrices (double S0, double T, double r, double sigma,double dt, int simSize);
    vector<vector<double>>
    genStockPricesForDelta(double S0, double T, double r, double sigma, double dt, double dS, int simsize);
    vector<double> genStockPricesT (double S0, double T, double r, double sigma,double dt, int simSize);
    double genStockPrices3 (double S0, double T, double r, double sigma, double dt, int simsize);
    void writeToFile (const string fNamePrefix, unique_ptr<double[]> v, const int size, const int k);
    void writeToFile (const string fNamePrefix, vector<double> x,vector<double> y,vector<double> ci, const int k);
    void writeToFile (const string fNamePrefix, vector<double> v, const int k, string colName);
};


#endif //MC_MONTECARLO_HPP
