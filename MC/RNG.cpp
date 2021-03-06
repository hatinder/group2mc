//
// Created by hsingh9 on 23/03/2019.
//

#include <iomanip>
#include <fstream>
#include <cassert>
#include <sstream>
#include "RNG.hpp"
#include "stats.hpp"

extern random_device rdStats;
using namespace stats;

//Box-Muller Method
vector<double> RNG::rngUsingBM (int size)
{
    vector<double> rng(size);
    mt19937 gen(rand());
    uniform_real_distribution<> distribution(0,1);
    for (int i = 0; i < size; i=i+2)
    {
        double r=distribution(gen);
        double phi=distribution(gen);
        rng[i]=sqrt(-2*log(r))*cos(2*M_PI*phi);
        rng[i+1]=sqrt(-2*log(r))*sin(2*M_PI*phi);
//        rng.push_back(sqrt(-2*log(r))*cos(2*M_PI*phi));
//        rng.push_back(sqrt(-2*log(r))*sin(2*M_PI*phi));
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

unique_ptr<double[]> RNG::rngUsingMTE (int size)
{
    unique_ptr<double[]> rng{new double[size]()};
    mt19937 gen(rand());
    uniform_real_distribution<> distribution(0.0,1.0);
    for (int i = 0; i < size; i=i+2)
    {
        double r=distribution(gen);
        double phi=distribution(gen);
        rng[i]=sqrt(-2*log(r))*cos(2*M_PI*phi);
        rng[i+1]=sqrt(-2*log(r))*sin(2*M_PI*phi);
    }
    return rng;
}

unique_ptr<double[]> RNG::rngUsingStatsBM (int size)
{
    unique_ptr<double[]> bm{new double[2]()};
//    random_device rd;
//    mt19937 gen(rand());
//    uniform_real_distribution<> distribution(0.0,1.0);
    double r=runif(0.0,1.0,rand());
    double phi=runif(0.0,1.0,rand());
//    double r=distribution(gen);
//    double phi=distribution(gen);
    bm[0]=sqrt(-2*log(r))*cos(2*M_PI*phi);
    bm[1]=sqrt(-2*log(r))*sin(2*M_PI*phi);
//    cout<<bm[0]<<endl;
//    cout<<bm[1]<<endl;
    return bm;
}

unique_ptr<double[]> RNG::rngUsingStatsBM2 (int size)
{
    unique_ptr<double[]> bm{new double[size]()};
    double r=runif(0.0,1.0,rand());
    double phi=runif(0.0,1.0,rand());
    for (int i = 0; i < size; i=i+2)
    {
        bm[i]=sqrt(-2*log(r))*cos(2*M_PI*phi);
        bm[i+1]=sqrt(-2*log(r))*sin(2*M_PI*phi);
    }
    return bm;
}

double RNG::rngUsingStatsBM ()
{
    unique_ptr<double[]> bm{new double[2]()};
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distribution(0.0,1.0);
//    double r=runif(0.0,1.0,rand());
//    double phi=runif(0.0,1.0,rand());
    double r=distribution(gen);
    double phi=distribution(gen);
    bm[0]=sqrt(-2*log(r))*cos(2*M_PI*phi);
    bm[1]=sqrt(-2*log(r))*sin(2*M_PI*phi);
//    cout<<bm[0]<<endl;
//    cout<<bm[1]<<endl;
    return bm[0];
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

void RNG::writeToFile (const string fNamePrefix, vector<double> v, const int k, string colName)
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
    for (int i = 0; i < v.size(); ++i)
    {
        oFileStream<<v[i]<<endl;
    }
    oFileStream.close();

}
