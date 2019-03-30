//
// Created by hsingh9 on 23/03/2019.
//

#include "RNG.hpp"

//Box-Muller Method
vector<double> RNG::rngUsingBM (int size)
{
    vector<double> rng(size);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distribution(0,1);
    for (int i = 0; i < size; i=i+2)
    {
        double r=distribution(gen);
        double phi=distribution(gen);
        rng[i]=sqrt(-2*log(r))*cos(2*M_PI*phi);
        rng[i+1]=sqrt(-2*log(r))*sin(2*M_PI*phi);

    }
//    for (int j = 0; j < rng.size(); ++j) {
//        cout<<rng[j]<<endl;
//    }
    return rng;
}

// Polar Marsaglia Method
vector<double> RNG::rngUsingPM ()
{
    vector<double> rng(2);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distribution(0,1);
    double w,u1,u2;
    do{
        u1=distribution(gen);
        u2=distribution(gen);
        w=u1*u1+u2*u2;
    }while(w>=1.0);
    w=sqrt((-2*log(w))/w);
    rng[0]=u1*w;
    rng[1]=u2*w;
    return rng;
}

vector<double> RNG::rngUsingMTE (int size)
{
    vector<double> rng(size);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distribution(0.0,1.0);
    for (int i = 0; i < rng.size(); ++i)
    {
        rng[i]=distribution(gen);
    }
    return rng;
}
