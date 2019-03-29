//
// Created by hsingh9 on 23/03/2019.
//

#ifndef MC_OPTIONDATA_HPP
#define MC_OPTIONDATA_HPP

#include <iostream>
#include <algorithm>
#include <memory>

using namespace std;


class OptionInfo
{
private:
    double S0, K, T, r, sigma, dy;
public:
    double getS0 () const;

    double getK () const;

    double getT () const;

    double getR () const;

    double getSigma () const;

    double getDY() const;
public:
    OptionInfo (double S0, double K, double T, double r, double sigma, double dy);
    OptionInfo();
    double payOff(double ST)
    {
        return max(ST-K,0.0);
    }

//    friend ostream& operator<< (ostream &out, const OptionInfo &oi);
    friend ostream& operator<< (ostream &out, const shared_ptr<OptionInfo> oi);

};

#endif //MC_OPTIONDATA_HPP
