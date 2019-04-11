//
// Created by Chris Broughall on 2019-03-28.
//


#include <cmath>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "EulerMethod.hpp"
#include "RNG.hpp"


vector<double> Euler::genStockPrices (double S0, double T, double r, double sigma,double dt, int simSize)
{
    //double dt = 0.5/simSize;
    int maxSize=1000000;
    double t;
    double a=dt*r, b=sigma * sqrt(dt);
    shared_ptr<RNG> rng=make_shared<RNG>();
    int tSize=round(T/dt);
    int rngSize=2*simSize*tSize;
//    rngSize = rngSize + (rngSize % 2 == 0 ? 2 : 1);
//    cout<<"simSize: "<<simSize<<" , tSize: "<<tSize <<", rngSize:"<< rngSize<<endl;
    if(rngSize<=maxSize)
    {
        vector<double> e = rng->rngUsingBM(rngSize);
//        cout<<"e size: "<<e.size()<<endl;
//        rng->writeToFile("SIMRNG",e,rngSize,"RNG");
        double St;
        vector<double> All(simSize, 0.0);
        int j = 0;
        for (int i = 0; i < simSize; i++)
        {
            t = 0.0;
            St=S0;
            do
            {
//                cout<<j<<endl;
                double val=e[j];
                St = St + (St * a) + (b * St * val);
//                cout<<"St:"<<St<<"e:"<<e[j]<<endl;
                t = t + dt;
                j++;
            } while (t <= T);
            All[i] = St;
        }
        return All;
    }
    else
    {
        int loopSize=floor(rngSize/maxSize);
        vector<int> sizeInfo(loopSize+1);
        int k=0;
        for (k = 0; k < loopSize; ++k)
        {
            sizeInfo[k]=maxSize;
        }
        sizeInfo[k]=rngSize-loopSize*maxSize;
        vector<double> e = rng->rngUsingBM(maxSize);
        double St = S0;
        vector<double> All(simSize, 0.0);
        int j = 0,si=1;
        for (int i = 0; i < simSize; i++)
        {
            t = 0.0;
            St = S0;
            do
            {

                St = St + (St * a) + (b * St * e[j]);
                t = t + dt;
                j++;
                if(j==maxSize)
                {
                    e=rng->rngUsingBM(sizeInfo[si++]);
                    j=0;
                }
            } while (t <= T);
            //std::cout << "The St is: " << St <<"\n";
            All[i] = St;
            }
            return All;
    }
}

vector<vector<double>> Euler::genStockPricesForDelta(double S0, double T, double r, double sigma, double dt, double dS,
                                                     int simSize)
{
    //double dt = 0.5/simSize;
    int maxsize=10000000;
    double t;
    double a=dt*r, b=sigma * sqrt(dt);
    shared_ptr<RNG> rng=make_shared<RNG>();
    int tSize=round(T/dt);
    int rngSize=simSize*tSize+simSize;
    if(rngSize<=maxsize)
    {
        vector<double> e = rng->rngUsingBM(rngSize);
//        cout<<"e size: "<<e.size()<<endl;
//        rng->writeToFile("SIMRNG",e,rngSize,"RNG");
        double StPlus, StMinus;
        vector<vector<double>> All;
        vector<double> AllPlus(simSize);
        vector<double> AllMinus(simSize);
        int j = 0;
        for (int i = 0; i < simSize; i++)
        {
            t = 0.0;
            StPlus = S0 + dS;
            StMinus = S0 - dS;
//            cout<<"StPlus: "<<StPlus<<endl;
//            cout<<"StMinus: "<<StMinus<<endl;
            do
            {
//                cout<<j<<endl;
                double val = e[j];
                StPlus = StPlus + (StPlus * a) + (b * StPlus * val);
                StMinus = StMinus + (StMinus * a) + (b * StMinus * val);
                t = t + dt;
                j++;
//                cout<<"St:"<<St<<"e:"<<e[j]<<endl;
            } while (t < T);
            AllPlus[i] = StPlus;
            AllMinus[i] = StMinus;
        }
        All.push_back(AllPlus);
        All.push_back(AllMinus);
        return All;
    }
    else
    {
        int loopSize=floor(rngSize/maxsize);
        vector<int> sizeInfo(loopSize+1);
        int k=0;
        for (; k < loopSize; ++k)
        {
            sizeInfo[k]=maxsize;
//            cout<<sizeInfo[k]<<endl;
        }
        //sizeInfo[k]=rngSize-loopSize*maxsize*1/2;
        sizeInfo[k]=maxsize;
        vector<double> e = rng->rngUsingBM(maxsize);
        double StPlus, StMinus;
        vector<vector<double>> All;
        vector<double> AllPlus(simSize);
        vector<double> AllMinus(simSize);
        int j = 0,si=1;
        for (int i = 0; i < simSize; i++)
        {
            t = 0.0;
            StPlus = S0 + dS;
            StMinus = S0 - dS;
            do {
//            St=St*exp(a+b*val);
                double val = e[j];
                StPlus = StPlus + (StPlus * a) + (b * StPlus * val);
                StMinus = StMinus + (StMinus * a) + (b * StMinus * val);
                t = t + dt;
                j++;
                if(j==maxsize)
                {
                    e.clear();
//                    cout<<"sizeinfo: "<<sizeInfo[si]<<endl;
                    e=rng->rngUsingBM(sizeInfo[si++]);
                    j=0;
                }
            } while (t < T);
//        cout<<"j: "<<j<<", "<<rngSize<<endl;
            AllPlus[i] = StPlus;
            AllMinus[i] = StMinus;
        }
        All.push_back(AllPlus);
        All.push_back(AllMinus);
        return All;

    }
}

double Euler::getStockPrice (double S0, double T, double r, double sigma, double dt)
{
    double St,StOld=S0;
    double t=0;
    double a=r*dt,b=sigma*sqrt(dt);
    unique_ptr<double[]> bmVal{new double[2]()};
    unique_ptr<RNG> rng;
    while (t<T)
    {
        bmVal= rng->rngUsingStatsBM(0);
        St=StOld + StOld*a + StOld*b*bmVal[0];
//        cout<<"St: "<<St<<endl;
        StOld=St;
        St=StOld + StOld*a + StOld*b*bmVal[1];
//        cout<<"St: "<<St<<endl;
        StOld=St;
        t=t+2*dt;
    }
//    bmVal.release();
    return St;
}


void Euler::writeToFile (const string fNamePrefix, unique_ptr<double[]> v, const int size, const int k)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    oFileStream<<setw(12)<<"x"<<setw(12)<<"y"<<endl;
    oFileStream<<endl;
    for (int i = 0; i < size; ++i)
    {
        oFileStream<<setw(12)<<i<<setw(12)<<v[i]<<endl;
    }
    oFileStream.close();
}


void Euler::writeToFile (const string fNamePrefix, unique_ptr<double[]> x,unique_ptr<double[]> y, const int size, const int k)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    oFileStream<<setw(12)<<"x"<<setw(12)<<"y"<<endl;
    oFileStream<<endl;
    for (int i = 0; i < size; ++i)
    {
        oFileStream<<setw(12)<<x[i]<<setw(12)<<y[i]<<endl;
    }
    oFileStream.close();
}

void Euler::writeToFile (const string fNamePrefix, unique_ptr<double[]> x,unique_ptr<double[]> y,unique_ptr<double[]> ci, const int size, const int k)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    oFileStream<<setw(12)<<"x"<<setw(12)<<"y"<<setw(12)<<"ci"<<endl;
    oFileStream<<endl;
    for (int i = 0; i < size; ++i)
    {
        oFileStream<<setw(12)<<x[i]<<setw(12)<<y[i]<<setw(12)<<ci[i]<<endl;
    }
    oFileStream.close();
}

void
Euler::writeToFile (const string fNamePrefix, vector<double> x, vector<double> y, vector<double> ci, const int k)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    oFileStream<<setw(12)<<"x"<<setw(12)<<"y"<<setw(12)<<"ci"<<endl;
    oFileStream<<endl;
    for (int i = 0; i < x.size(); ++i)
    {
        oFileStream<<setw(12)<<x[i]<<setw(12)<<y[i]<<setw(12)<<ci[i]<<endl;
    }
    oFileStream.close();

}
