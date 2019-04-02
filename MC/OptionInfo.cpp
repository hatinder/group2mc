//
// Created by hsingh9 on 23/03/2019.
//

#include "OptionInfo.hpp"



OptionInfo::OptionInfo (double S0, double K, double T, double r, double sigma) : S0(S0), K(K), T(T), r(r), sigma(sigma)
{
}

OptionInfo::OptionInfo (): S0(70.0), K(65.0), T(0.5), r(0.07), sigma(0.27)
{

}

double OptionInfo::getS0 () const
{
    return S0;
}

double OptionInfo::getK () const
{
    return K;
}

double OptionInfo::getT () const
{
    return T;
}

double OptionInfo::getR () const
{
    return r;
}

double OptionInfo::getSigma () const
{
    return sigma;
}

//ostream& operator<< (std::ostream &out, const OptionInfo &oi)
//{
////    out <<endl<<"Option( S0: " << oi.S0 << ", K" << oi.K << ", T" << oi.T << "r: " << oi.r <<"sigma: " << oi.sigma<< ")"<<endl;
//
//    return out;
//}
//
ostream& operator<< (std::ostream &out, const shared_ptr<OptionInfo> oi)
{
    out <<"Option( S0: " << oi->S0 << ", K: " << oi->K << ", T: " << oi->T << ", r: " << oi->r <<", sigma: " << oi->sigma<< ")";
    return out;
}
