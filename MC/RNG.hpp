//
// Created by hsingh9 on 23/03/2019.
//

#ifndef MC_RANDOMNUMBERGENERATOR_HPP
#define MC_RANDOMNUMBERGENERATOR_HPP

#include <vector>
#include <iostream>
#include <random>
#include <memory>

using namespace std;

class RNG
{
public:
    vector<double> rngUsingBM (int size);    //Box Muller Method
    vector<double> rngUsingPM();    //Polar Marsaglia Method
    vector<double> rngUsingMTE (int size);   //Mersenne Twister Engine
    unique_ptr<double[]> rngUsingStatsBM (int size);
    void writeToFile (const string fNamePrefix, unique_ptr<double[]> v, const int size, const int k, string val);

};


#endif //MC_RANDOMNUMBERGENERATOR_HPP
