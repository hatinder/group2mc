//
// Created by hsingh9 on 23/03/2019.
//

#include <iomanip>
#include <fstream>
#include <cassert>
#include <sstream>
#include "RNG.hpp"
#include "stats.hpp"

using namespace stats;

//Box-Muller Method
vector<double> RNG::rngUsingBM (int size)
{
    vector<double> rng(size);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distribution(0,1);
    for (int i = 0; i < size-1; i=i+2)
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

unique_ptr<double[]> RNG::rngUsingStatsBM (int size)
{
    unique_ptr<double[]> bm{new double[2]()};
//    random_device rd;
    mt19937 gen(rand());
    uniform_real_distribution<> distribution(0.0,1.0);
//    double val=runif(0.0,1.0);
    double val=distribution(gen);
    bm[0]=sqrt(-2*log(val))*cos(2*M_PI*val);
    bm[1]=sqrt(-2*log(val))*sin(2*M_PI*val);
    cout<<bm[0]<<endl;
//    cout<<bm[1]<<endl;
    return bm;
}

void RNG::writeToFile (const string fNamePrefix, unique_ptr<double[]> v, const int size, const int k, string colName)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    oFileStream<<setw(12)<<colName;
    oFileStream<<endl;
    for (int i = 0; i < size; ++i)
    {
        oFileStream<<v[i]<<endl;
    }
    oFileStream.close();

}
