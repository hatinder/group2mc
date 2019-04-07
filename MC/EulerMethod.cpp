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


vector<double> Euler::genStockPrices (double S0, double T, double r, double sigma, int simSize)
{

    //double dt = 0.5/simSize;
    double dt = 1.0/365.0;
    double t;



    double St = S0;
    vector<double> All(simSize,0.0);
    double StOld = St;

    int j = 0;

        for (int i = 0; i < simSize; i++) {

            t = 0.0;
            St = S0;
            shared_ptr<RNG> rng=make_shared<RNG>();
            vector<double> e= rng->rngUsingBM(simSize*round(T/dt) + 1);
            j = 0;
            StOld = S0;

            do{

                St = StOld + (StOld * dt*r) + (sigma * sqrt(dt) * StOld * e[j]);
                StOld = St;


                t = t + dt;
                j++;

            }while(t <= T);

            //std::cout << "The St is: " << St <<"\n";

            All[i] = St;

        }



    return All;

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
