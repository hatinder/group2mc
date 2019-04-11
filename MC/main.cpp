#include <iostream>
#include <memory>
#include "OptionInfo.hpp"
#include "RNG.hpp"
#include "MonteCarlo.hpp"
#include "EulerMethod.hpp"
#include "BlackScholes.hpp"
#include <iomanip>
#include <fstream>


random_device rdStats;

void initOptionInfo();
void runMonteCarloSimulation (shared_ptr<OptionInfo> optionInfo, int MCSIM);
void runRNG (int size);
void runEuler (shared_ptr<OptionInfo> optionInfo, int simSize);
void runBlackScholes();
void runEulerAssetOrNothing();
void runMonteCarloSimulationForAssetOrNothing (shared_ptr<OptionInfo> optionInfo, int MCSIM);
void runAssetOrNothing();
void runRNGUsingStatsBM();
void runRNGUsingStatsBM1 (int size);
void runWeakEulerStockPrice();
void runWeakEulerMCCP();
void runMCUsingStatsBM();
void runWeakEulerConvergenceWithDt();
void runWeakEulerConvergenceMCError();
void runRNGUsingStatsBM2();
void runMCDelta (shared_ptr<OptionInfo> optionInfo, int sSize);
void runMCAONDelta (shared_ptr<OptionInfo> optionInfo, int MCSIM);
void runEulerDelta (shared_ptr<OptionInfo> optionInfo, int simSize);
void runEulerAON (shared_ptr<OptionInfo> optionInfo, int simSize);
void runEulerDeltaAON (shared_ptr<OptionInfo> optionInfo, int simSize);
void runEulerDeltaFDM (shared_ptr<OptionInfo> optionInfo, int simSize);
void runEulerDeltaAONFDM (shared_ptr<OptionInfo> optionInfo, int simSize);
bool run();

int main (int argc, char *argv[])
{
    bool val=run();
    if(val)
        cout<<"Completed Successfully";
    else
        cout<<"Run time problem. Please check";
//    int sSize=1000;
//    if(argc==2)
//    {
//        sSize=atoi(argv[1]);
//        cout<<"Starting with simSize: "<<sSize<<endl;
//    }
//    initOptionInfo();
//    runBlackScholes();

//    runMonteCarloSimulation();
//    runMCDelta(sSize);
//        runEuler();
//    cout<<endl;
//    runAssetOrNothing();
//    runEulerAON();
//      runEulerDeltaAON();
//     runEulerDeltaFDM();
//    runEulerDelta();
//      runEulerDeltaAONFDM();
//    runMonteCarloSimulationForAssetOrNothing();
//    runMCAONDelta(sSize);
//    runEulerAssetOrNothing();
//    runRNGUsingStatsBM();
//    clock_t begin=clock();
//    runRNGUsingStatsBM1(10000000);
//    clock_t end=clock();
//    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    cout<<"Time Taken UP: "<<elapsed_secs<<endl;
//    runWeakEulerStockPrice();
//    runWeakEulerMCCP();
//    runMCUsingStatsBM();
//    begin=clock();
//    runRNG(1000);
//    end=clock();
//    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    cout<<"Time Taken vector: "<<elapsed_secs<<endl;
//    runWeakEulerConvergenceWithDt();
//    runWeakEulerConvergenceMCError();
//    runRNGUsingStatsBM2();
    return 0;
}

bool run()
{
    ifstream ifs;
    string filename="config.txt";
    ifs.open(filename.c_str(), ios_base::in);
    if (!ifs)
    {
        cout << "config.txt file is missing: " << filename << endl;
        return false;
    }
    vector<double> parameterValues;
    vector<string> parameterNames ;
    double value;
    string name;
//    ifs >> name;
//    cout << name<<endl;
    while (ifs >> name >> value)
    {
        parameterNames.push_back(name);
        parameterValues.push_back(value);
    }
//    cout<<"config file parameters"<<endl;
//    for (int i = 0; i < parameterValues.size(); ++i)
//    {
//        cout<<"Name: "<<parameterNames[i]<<" , Value: "<<parameterValues[i]<<endl;
//    }
    ifs.close();
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3], parameterValues[4]);
    cout<<  optionInfo<<endl;
    int sampleSize=parameterValues[5];
    int simSize=parameterValues[6];
    cout<<"Sample Size: "<<sampleSize<<", Simulation Size: "<<simSize<<endl;

    cout<<"=========== "<< endl;
    cout<<"1.Running European Call for Monte Carlo Error with Sample Size: "<<sampleSize<< endl;
    cout<<"=========== "<< endl;
    runMonteCarloSimulation(optionInfo,sampleSize);

    cout<<"=========== "<< endl;
    cout<<"2.Running European Call for Euler Error with simulation Size: "<<simSize<< endl;
    cout<<"=========== "<< endl;
    runEuler(optionInfo,simSize);

    cout<<"=========== "<< endl;
    cout<<"3.Running European Call Delta for Monte Carlo Error with Sample Size: "<<sampleSize<< endl;
    cout<<"=========== "<< endl;
    runMCDelta(optionInfo,sampleSize);

    cout<<"=========== "<< endl;
    cout<<"4.Running European Call Delta for Euler Error with simulation Size: "<<simSize<< endl;
    cout<<"=========== "<< endl;
    runEulerDelta(optionInfo,simSize);

    cout<<"=========== "<< endl;
    cout<<"5.Running European Call Delta for Finite Difference Method (dS) error with simulation Size: "<<simSize<< endl;
    cout<<"=========== "<< endl;
    runEulerDeltaFDM(optionInfo,simSize);

    cout<<"=========== "<< endl;
    cout<<"6.Running Asset-Or-Nothing Call Monte Carlo Error sample size: "<<sampleSize<< endl;
    cout<<"=========== "<< endl;
    runMonteCarloSimulationForAssetOrNothing (optionInfo,sampleSize);

    cout<<"=========== "<< endl;
    cout<<"7.Running Asset-Or-Nothing Call Euler Error with simulation Size: "<<simSize<< endl;
    cout<<"=========== "<< endl;
    runEulerAON(optionInfo,simSize);

    cout<<"=========== "<< endl;
    cout<<"8.Running Asset-Or-Nothing Call Delta Monte Carlo Error with Sample Size: "<<sampleSize<< endl;
    cout<<"=========== "<< endl;
    runMCAONDelta(optionInfo,sampleSize);

    cout<<"=========== "<< endl;
    cout<<"9.Running Asset-Or-Nothing Call Delta Euler Error with simulation Size: "<<simSize<< endl;
    cout<<"=========== "<< endl;
    runEulerDeltaAON(optionInfo,simSize);

    cout<<"=========== "<< endl;
    cout<<"10.Running Asset-Or-Nothing Call Delta for Finite Difference Method (dS) error with simulation Size: "<<simSize<< endl;
    cout<<"=========== "<< endl;
    runEulerDeltaAONFDM(optionInfo,simSize);

    return true;
}

void initOptionInfo()
{
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    cout<<  optionInfo<<endl;
//    shared_ptr<OptionInfo> optionInfo2=make_shared<OptionInfo>(70.0, 65.0, 0.5, 0.07, 0.27);
//    cout<<  optionInfo2<<endl;

}

void runEuler (shared_ptr<OptionInfo> optionInfo, int simSize)
{
    srand(2019);
    unique_ptr<Euler> euler;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int dtSize = 6;
    double dt[dtSize] = {0.5, 0.25, 0.125, 0.0625,0.015625,0.0078125};
    vector<double> x(dtSize),y(dtSize),ci(dtSize);
//    int simSize = 10000000;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());
//    ,SjMinusDT(simSize,0.0),SjPlusDT(simSize,0.0);
//    vector<double> VjMinusDT(simSize,0.0),VjPlusDT(simSize,0.0);
    int p=0;
    for (int i = 0; i < dtSize; ++i)
    {
        vector<double> Sj(simSize), Vj(simSize);
        double sigmaHatSqr=0.0;
        clock_t begin = clock();
        Sj = euler->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(),
                                   dt[i], simSize);
//        cout << "SJ Size:" << Sj.size() << endl;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        for (int j = 0; j < Sj.size(); ++j)
        {
            Vj[j] = optionInfo->payOff(Sj[j]);
        }
//        double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//        cout << "Approximate Stock Price at T : " << avgST << endl;
        double avgVT = accumulate(Vj.begin(), Vj.end(), 0.0) / Vj.size();
        double eulerCallPrice=D * avgVT;
        for (int k = 0; k < Vj.size(); ++k)
        {
            sigmaHatSqr += (avgVT - Vj[k]) * (avgVT - Vj[k]);
        }
        double avgSigmaHat =sigmaHatSqr/Vj.size();
        cout << "Euler Call Price: " << eulerCallPrice  << " with CI: (+-) " << D * 1.96 * sqrt(avgSigmaHat / Vj.size()) << ", dt: " <<dt[i]<<endl;
        x[i]=dt[i];
        y[i]=values[0]-eulerCallPrice;
        ci[i]=D * 1.96 * sqrt(avgSigmaHat / Vj.size());
    }
    euler->writeToFile("EULER",x,y,ci,1);
}

void runEulerDelta (shared_ptr<OptionInfo> optionInfo, int simSize)
{
    srand(2019);
    unique_ptr<Euler> euler;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int dtSize = 6;
    double dt[dtSize] = {0.5, 0.25, 0.125, 0.0625,0.015625,0.0078125};
    vector<double> x(dtSize),y(dtSize),ci(dtSize);
//    int simSize = 1000000;
    double t = 0.0,dS=1.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());

//    cout<<"BS: "<<values[0]<<endl;
     double avgMCPricePlus=0.0,avgMCPriceMinus=0.0;
    clock_t begin = clock();
    double MCPlus,MCMinus,avgSigmaHat;
    vector<double> eulerDelta(dtSize,0.0);
    for (int l = 0; l < dtSize; ++l) {
        vector<vector<double>> Sj;
        double payOffMinus = 0.0, payOffPlus = 0.0;
        Sj = euler->genStockPricesForDelta(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),
                                           optionInfo->getSigma(), dt[l], dS, simSize);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
//            cout<<"Sj[0][0]"<<Sj[0][0]<<endl;
//            cout<<"Sj[1][0]"<<Sj[1][0]<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffPlus += optionInfo->payOff(Sj[0][j]);

        }
//            cout<<"Pay Off Plus: "<<payOffPlus<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffMinus += optionInfo->payOff(Sj[1][j]);

        }
//            cout<<"Pay Off Minus: "<<payOffMinus<<endl;
        double avgVTPlus = payOffPlus / simSize;
        double avgVTMinus = payOffMinus / simSize;
//        cout<<"avgVTPlus: "<<avgVTPlus<<endl;
//        cout<<"avgVTMinus: "<<avgVTMinus<<endl;
        MCPlus = D * avgVTPlus;
        MCMinus = D * avgVTMinus;
//        cout <<"MCPlus: "<< MCPlus<<endl;
//        cout <<"MCMinus: "<<MCMinus<<endl;
        eulerDelta[l] = (MCPlus - MCMinus) / (2.0 * dS);

//            cout<<"avgPricePlus"<<avgMCPricePlus<<endl;
//            cout<<"avgPriceMinus"<<avgMCPriceMinus<<endl;
        double sigmaHatSqr = 0.0;
        for (int k = 0; k < simSize; ++k) {
            double eulerKDelta=(optionInfo->payOff(Sj[0][k])-optionInfo->payOff(Sj[1][k]))/(2.0*dS);
            sigmaHatSqr += (eulerDelta[l] - eulerKDelta) * (eulerDelta[l] - eulerKDelta);
        }
        avgSigmaHat = sigmaHatSqr / simSize;
        ci[l] += D * 1.96 * sqrt(avgSigmaHat / simSize);
        y[l]=values[1]-eulerDelta[l];
        x[l]=dt[l];
//        cout<<"y before: "<<y[i]<<endl;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        cout << "MC Delta : " << eulerDelta[l]<< " with CI: (+-) " << ci[l] << endl;
    }
    euler->writeToFile("EULERDELTA",x,y,ci,1);
}

void runEulerDeltaFDM (shared_ptr<OptionInfo> optionInfo, int simSize)
{
    srand(2019);
    unique_ptr<Euler> euler;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    double dt = 1.0/365.0;
    const int dSSize = 3;
    double dS[dSSize] = {5.0,4.0,3.0};
    vector<double> x(dSSize),y(dSSize),ci(dSSize);
//    int simSize = 1000000;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());

//    cout<<"BS: "<<values[0]<<endl;
    double avgMCPricePlus=0.0,avgMCPriceMinus=0.0;
    clock_t begin = clock();
    double MCPlus,MCMinus,avgSigmaHat;
    vector<double> eulerDelta(dSSize,0.0);
    for (int l = 0; l < dSSize; ++l) {
        vector<vector<double>> Sj;
        double payOffMinus = 0.0, payOffPlus = 0.0;
        Sj = euler->genStockPricesForDelta(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),
                                           optionInfo->getSigma(), dt, dS[l], simSize);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
//            cout<<"Sj[0][0]"<<Sj[0][0]<<endl;
//            cout<<"Sj[1][0]"<<Sj[1][0]<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffPlus += optionInfo->payOff(Sj[0][j]);

        }
//            cout<<"Pay Off Plus: "<<payOffPlus<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffMinus += optionInfo->payOff(Sj[1][j]);

        }
//            cout<<"Pay Off Minus: "<<payOffMinus<<endl;
        double avgVTPlus = payOffPlus / simSize;
        double avgVTMinus = payOffMinus / simSize;
//        cout<<"avgVTPlus: "<<avgVTPlus<<endl;
//        cout<<"avgVTMinus: "<<avgVTMinus<<endl;
        MCPlus = D * avgVTPlus;
        MCMinus = D * avgVTMinus;
//        cout <<"MCPlus: "<< MCPlus<<endl;
//        cout <<"MCMinus: "<<MCMinus<<endl;
        eulerDelta[l] = (MCPlus - MCMinus) / (2.0 * dS[l]);

//            cout<<"avgPricePlus"<<avgMCPricePlus<<endl;
//            cout<<"avgPriceMinus"<<avgMCPriceMinus<<endl;
        double sigmaHatSqr = 0.0;
        for (int k = 0; k < simSize; ++k) {
            double eulerKDelta=(optionInfo->payOff(Sj[0][k])-optionInfo->payOff(Sj[1][k]))/(2.0*dS[l]);
            sigmaHatSqr += (eulerDelta[l] - eulerKDelta) * (eulerDelta[l] - eulerKDelta);
        }
        avgSigmaHat = sigmaHatSqr / simSize;
        ci[l] += D * 1.96 * sqrt(avgSigmaHat / simSize);
        y[l]=values[1]-eulerDelta[l];
        x[l]=dS[l];
//        cout<<"y before: "<<y[i]<<endl;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        cout << "MC Delta : " << eulerDelta[l]<< ", dS: " << dS[l]
             << " with CI: (+-) " << ci[l] << endl;
    }
    euler->writeToFile("EULERDELTAFDM",x,y,ci,1);
}

void runEulerDeltaAONFDM (shared_ptr<OptionInfo> optionInfo, int simSize)
{
    srand(2019);
    unique_ptr<Euler> euler;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    double dt = 1.0/365.0;
    const int dSSize = 3;
    double dS[dSSize] = {3.0,2.0,1.0};
    vector<double> x(dSSize),y(dSSize),ci(dSSize);
//    int simSize = 1000000;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values= bs->calcAssetOrNothingExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());

//    cout<<"BS: "<<values[0]<<endl;
    double avgMCPricePlus=0.0,avgMCPriceMinus=0.0;
    clock_t begin = clock();
    double MCPlus,MCMinus,avgSigmaHat;
    vector<double> eulerDelta(dSSize,0.0);
    for (int l = 0; l < dSSize; ++l) {
        vector<vector<double>> Sj;
        double payOffMinus = 0.0, payOffPlus = 0.0;
        Sj = euler->genStockPricesForDelta(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),
                                           optionInfo->getSigma(), dt, dS[l], simSize);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
//            cout<<"Sj[0][0]"<<Sj[0][0]<<endl;
//            cout<<"Sj[1][0]"<<Sj[1][0]<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffPlus += optionInfo->payOffForAssetOrNothing(Sj[0][j]);

        }
//            cout<<"Pay Off Plus: "<<payOffPlus<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffMinus += optionInfo->payOffForAssetOrNothing(Sj[1][j]);

        }
//            cout<<"Pay Off Minus: "<<payOffMinus<<endl;
        double avgVTPlus = payOffPlus / simSize;
        double avgVTMinus = payOffMinus / simSize;
//        cout<<"avgVTPlus: "<<avgVTPlus<<endl;
//        cout<<"avgVTMinus: "<<avgVTMinus<<endl;
        MCPlus = D * avgVTPlus;
        MCMinus = D * avgVTMinus;
//        cout <<"MCPlus: "<< MCPlus<<endl;
//        cout <<"MCMinus: "<<MCMinus<<endl;
        eulerDelta[l] = (MCPlus - MCMinus) / (2.0 * dS[l]);

//            cout<<"avgPricePlus"<<avgMCPricePlus<<endl;
//            cout<<"avgPriceMinus"<<avgMCPriceMinus<<endl;
        double sigmaHatSqr = 0.0;
        for (int k = 0; k < simSize; ++k) {
            double eulerKDelta=(optionInfo->payOffForAssetOrNothing(Sj[0][k])-optionInfo->payOffForAssetOrNothing(Sj[1][k]))/(2.0*dS[l]);
            sigmaHatSqr += (eulerDelta[l] - eulerKDelta) * (eulerDelta[l] - eulerKDelta);
        }
        avgSigmaHat = sigmaHatSqr / simSize;
        ci[l] += D * 1.96 * sqrt(avgSigmaHat / simSize);
        y[l]=values[1]-eulerDelta[l];
        x[l]=dS[l];
//        cout<<"y before: "<<y[i]<<endl;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        cout << "MC Delta : " << eulerDelta[l]<< ", dS: " << dS[l]
             << " with CI: (+-) " << ci[l] << endl;
    }
    euler->writeToFile("EULERDELTAAONFDM",x,y,ci,1);
}

void runEulerAON (shared_ptr<OptionInfo> optionInfo, int simSize)
{
    srand(2019);
    unique_ptr<Euler> euler;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int dtSize = 6;
    double dt[dtSize] = {0.5, 0.25, 0.125, 0.0625,0.015625,0.0078125};
    vector<double> x(dtSize),y(dtSize),ci(dtSize);
//    int simSize = 10000000;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values= bs->calcAssetOrNothingExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());
//    ,SjMinusDT(simSize,0.0),SjPlusDT(simSize,0.0);
//    vector<double> VjMinusDT(simSize,0.0),VjPlusDT(simSize,0.0);
    int p=0;
    for (int i = 0; i < dtSize; ++i)
    {
        vector<double> Sj(simSize), Vj(simSize);
        double sigmaHatSqr=0.0;
        clock_t begin = clock();
        Sj = euler->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(),
                                   dt[i], simSize);
//        cout << "SJ Size:" << Sj.size() << endl;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        for (int j = 0; j < Sj.size(); ++j)
        {
            Vj[j] = optionInfo->payOffForAssetOrNothing(Sj[j]);
        }
//        double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//        cout << "Approximate Stock Price at T : " << avgST << endl;
        double avgVT = accumulate(Vj.begin(), Vj.end(), 0.0) / Vj.size();
        double eulerCallPrice=D * avgVT;
        for (int k = 0; k < Vj.size(); ++k)
        {
            sigmaHatSqr += (avgVT - Vj[k]) * (avgVT - Vj[k]);
        }
        double avgSigmaHat =sigmaHatSqr/Vj.size();
        cout << "Euler Call Price: " << eulerCallPrice  << " with CI: (+-) " << D * 1.96 * sqrt(avgSigmaHat / Vj.size()) << " for dt: "<<dt[i]<< endl;
        x[i]=dt[i];
        y[i]=values[0]-eulerCallPrice;
        ci[i]=D * 1.96 * sqrt(avgSigmaHat / Vj.size());
    }
    euler->writeToFile("EULERAON",x,y,ci,1);
}

void runEulerDeltaAON (shared_ptr<OptionInfo> optionInfo, int simSize)
{
    srand(2019);
    unique_ptr<Euler> euler;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int dtSize = 5;
    double dt[dtSize] = {0.5, 0.25, 0.125, 0.0625,0.015625};
    vector<double> x(dtSize),y(dtSize),ci(dtSize);
//    int simSize = 10000000;
    double t = 0.0,dS=1.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values= bs->calcAssetOrNothingExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());

//    cout<<"BS AON CP: "<<values[0]<<endl;
//    cout<<"BS AON DELTA: "<<values[1]<<endl;
    double avgMCPricePlus=0.0,avgMCPriceMinus=0.0;
    clock_t begin = clock();
    double MCPlus,MCMinus,avgSigmaHat;
    vector<double> eulerDelta(dtSize,0.0);
    for (int l = 0; l < dtSize; ++l) {
        vector<vector<double>> Sj;
        double payOffMinus = 0.0, payOffPlus = 0.0;
        Sj = euler->genStockPricesForDelta(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),
                                           optionInfo->getSigma(), dt[l], dS, simSize);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
//            cout<<"Sj[0][0]"<<Sj[0][0]<<endl;
//            cout<<"Sj[1][0]"<<Sj[1][0]<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffPlus += optionInfo->payOffForAssetOrNothing(Sj[0][j]);

        }
//            cout<<"Pay Off Plus: "<<payOffPlus<<endl;
        for (int j = 0; j < simSize; ++j) {
            payOffMinus += optionInfo->payOffForAssetOrNothing(Sj[1][j]);

        }
//            cout<<"Pay Off Minus: "<<payOffMinus<<endl;
        double avgVTPlus = payOffPlus / simSize;
        double avgVTMinus = payOffMinus / simSize;
//        cout<<"avgVTPlus: "<<avgVTPlus<<endl;
//        cout<<"avgVTMinus: "<<avgVTMinus<<endl;
        MCPlus = D * avgVTPlus;
        MCMinus = D * avgVTMinus;
//        cout <<"MCPlus: "<< MCPlus<<endl;
//        cout <<"MCMinus: "<<MCMinus<<endl;
        eulerDelta[l] = (MCPlus - MCMinus) / (2.0 * dS);

//            cout<<"avgPricePlus"<<avgMCPricePlus<<endl;
//            cout<<"avgPriceMinus"<<avgMCPriceMinus<<endl;
        double sigmaHatSqr = 0.0;
        for (int k = 0; k < simSize; ++k) {
            double eulerKDelta=(optionInfo->payOffForAssetOrNothing(Sj[0][k])-optionInfo->payOffForAssetOrNothing(Sj[1][k]))/(2.0*dS);
            sigmaHatSqr += (eulerDelta[l] - eulerKDelta) * (eulerDelta[l] - eulerKDelta);
        }
        avgSigmaHat = sigmaHatSqr / simSize;
        ci[l] += D * 1.96 * sqrt(avgSigmaHat / simSize);
        y[l]=values[1]-eulerDelta[l];
        x[l]=dt[l];
//        cout<<"y before: "<<y[i]<<endl;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        cout << "MC Delta : " << eulerDelta[l]<< ", simSize: " << simSize
             << " with CI: (+-) " << ci[l] << endl;
    }
    euler->writeToFile("EULERDELTAAON",x,y,ci,1);
}

void runMonteCarloSimulation (shared_ptr<OptionInfo> optionInfo, int MCSIM)
{
    srand(2019);
    unique_ptr<MonteCarlo> mc;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int NSIM = 5;
//    const int MCSIM=1000;
    double simSize[NSIM] = {100,200,300,400,500};
    vector<double> x(NSIM,0.0),y(NSIM,0.0),ci(NSIM,0.0);
    double dt = 1.0 / 365.0;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values = bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(),
                                                optionInfo->getSigma(), optionInfo->getT());
//    cout<<"BS: "<<values[0]<<endl;
    for (int i = 0; i < NSIM; ++i)
    {
        double avgMCPrice=0.0;
        clock_t begin = clock();
        double MCPrice,avgSigmaHat;
        vector<double> mcPrice(MCSIM,0.0),mcError(MCSIM,0.0);
        for (int l = 0; l < MCSIM; ++l)
        {
            vector<double> Sj(simSize[i], 0.0);
            double payOff = 0.0;
            Sj = mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(), dt,
                                    simSize[i]);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
            for (int j = 0; j < Sj.size(); ++j)
            {
                payOff += optionInfo->payOff(Sj[j]);

            }
            double avgVT = payOff / Sj.size();
            double sigmaHatSqr = 0.0;
            for (int k = 0; k < Sj.size(); ++k)
            {
                sigmaHatSqr += (avgVT - optionInfo->payOff(Sj[k])) * (avgVT - optionInfo->payOff(Sj[k]));
            }
            avgSigmaHat = sigmaHatSqr / Sj.size();
            MCPrice=D*avgVT;
            mcPrice[l]=MCPrice;
            mcError[l]=(values[0]-MCPrice);
//            cout << MCPrice<<endl;
            avgMCPrice += MCPrice;
            y[i]+=(values[0]-MCPrice);
            ci[i]+=D * 1.96 * sqrt(avgSigmaHat / Sj.size());
        }
        mc->writeToFile("MCPRICE",mcPrice,i,"MCP");
        mc->writeToFile("MCERR",mcError,i,"ERR");
//        cout<<"y before: "<<y[i]<<endl;
        x[i]=sqrt(1.0/simSize[i]);
        y[i]=y[i]/MCSIM;
        ci[i]=ci[i]/MCSIM;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize[i] << " is :" << elapsed_secs << endl;
        cout << "MC Call Price: " << avgMCPrice/MCSIM <<", simSize: "<< simSize[i] <<" with CI: (+-) " << ci[i] << endl;

    }
    mc->writeToFile("MCERROR",x,y,ci,1);
}


void runMCDelta (shared_ptr<OptionInfo> optionInfo, int MCSIM)
{
    srand(2019);
    unique_ptr<MonteCarlo> mc;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int NSIM = 4;
//    const int MCSIM = 1000;
    int sSize=100;
    int simSize;
    vector<double> x(NSIM,0.0),y(NSIM,0.0),ci(NSIM,0.0);
    double dt = 1.0 / 365.0;
    double dS=1.0;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values = bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(),
                                                optionInfo->getSigma(), optionInfo->getT());
//    cout<<"BS: "<<values[0]<<endl;
    for (int i = 0; i < NSIM; ++i)
    {
        simSize=(i+1)*sSize;
        double avgMCPricePlus=0.0,avgMCPriceMinus=0.0;
        clock_t begin = clock();
        double MCPlus,MCMinus,avgSigmaHat;
        vector<double> mcDelta(MCSIM,0.0),mcError(MCSIM,0.0);
        for (int l = 0; l < MCSIM; ++l)
        {
            vector<vector<double>> Sj;
            double payOffMinus = 0.0,payOffPlus=0.0;
            Sj = mc->genStockPricesForDelta(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),
                                            optionInfo->getSigma(), dt, dS,
                                            simSize);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
//            cout<<"Sj[0][0]"<<Sj[0][0]<<endl;
//            cout<<"Sj[1][0]"<<Sj[1][0]<<endl;
            for (int j = 0; j < simSize; ++j)
            {
                payOffPlus += optionInfo->payOff(Sj[0][j]);

            }
//            cout<<"Pay Off Plus: "<<payOffPlus<<endl;
            for (int j = 0; j < simSize; ++j)
            {
                payOffMinus += optionInfo->payOff(Sj[1][j]);

            }
//            cout<<"Pay Off Minus: "<<payOffMinus<<endl;
            double avgVTPlus = payOffPlus / simSize;
            double avgVTMinus = payOffMinus / simSize;
//            cout<<"avgVTPlus: "<<avgVTPlus<<endl;
//            cout<<"avgVTMinus: "<<avgVTMinus<<endl;
            MCPlus=D*avgVTPlus;
            MCMinus=D*avgVTMinus;
//            cout << MCPrice<<endl;
            avgMCPricePlus += MCPlus;
            avgMCPriceMinus += MCMinus;
            mcDelta[l]=(MCPlus-MCMinus)/(2*dS);
            mcError[l]=values[1]-mcDelta[l];
//            cout<<"avgPricePlus"<<avgMCPricePlus<<endl;
//            cout<<"avgPriceMinus"<<avgMCPriceMinus<<endl;
            double sigmaHatSqr = 0.0;
            for (int k = 0; k < simSize; ++k)
            {
                sigmaHatSqr += (avgVTPlus - optionInfo->payOff(Sj[0][k])) * (avgVTPlus - optionInfo->payOff(Sj[0][k]));
            }
            avgSigmaHat = sigmaHatSqr / Sj.size();
            ci[i]+=D * 1.96 * sqrt(avgSigmaHat / simSize);
        }
//        cout<<"y before: "<<y[i]<<endl;
        mc->writeToFile("MCDELTA",mcDelta,i,"MCP");
        mc->writeToFile("MCDELTAERR",mcError,i,"ERR");
        x[i]=sqrt(1.0/simSize);
        y[i]=y[i]/MCSIM;
        double avgMCDelta = accumulate(mcDelta.begin(), mcDelta.end(), 0.0) / mcDelta.size();
        double tempCI=0.0;
        for (int m = 0; m < MCSIM; ++m)
        {
            tempCI += (avgMCDelta-mcDelta[m])*(avgMCDelta-mcDelta[m]);
        }
        ci[i]=D*1.96*tempCI/MCSIM;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        cout << "MC Delta : " << (avgMCPricePlus-avgMCPriceMinus)/(MCSIM*2*dS) <<", simSize: "<< simSize <<" with CI: (+-) " << ci[i] << endl;

    }
    mc->writeToFile("MCERROR",x,y,ci,1);
}

    //DELTA SIMULATION
//    SjMinusDT= mc->genStockPrices3(optionInfo->getS0(), optionInfo->getT() - dt, optionInfo->getR(),
//                                   optionInfo->getSigma(), simSize);
//    SjPlusDT = mc->genStockPrices3(optionInfo->getS0(), optionInfo->getT() + dt, optionInfo->getR(),
//                                   optionInfo->getSigma(), simSize);
//    for (int j = 0; j < Sj.simSize(); ++j)
//    {
//                cout<<Sj[j]<<endl;
//        VjMinusDT[j]=optionInfo->payOff(SjMinusDT[j]);
//        VjPlusDT[j]=optionInfo->payOff(SjPlusDT[j]);
//    }
//
//    double avgVTMinusDT = accumulate( VjMinusDT.begin(), VjMinusDT.end(), 0.0)/VjMinusDT.simSize();
//    double avgVTPlusDT = accumulate( VjPlusDT.begin(), VjPlusDT.end(), 0.0)/VjPlusDT.simSize();
//    cout << "MC Delta: " << ((exp(-optionInfo->getR()*(optionInfo->getT()+dt-t))*avgVTPlusDT)-(exp(-optionInfo->getR()*(optionInfo->getT()-dt-t))*avgVTMinusDT))/(2*dt) << endl;

//}

void runBlackScholes()
{
    unique_ptr<BlackScholes> bs;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                  optionInfo->getT());
    cout<<"BS Call Price: "<<values[0]<<endl;
    cout<<"BS Delta: "<<values[1]<<endl;

}

void runEulerAssetOrNothing()
{
    unique_ptr<Euler> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=100;
    double t = 0.0,dt=1.0/365.0;

    vector<double> Sj(size,0.0),Vj(size,0.0);
    Sj = mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),dt,size);
    for (int j = 0; j < Sj.size(); ++j)
    {
        Vj[j]=optionInfo->payOffForAssetOrNothing(Sj[j]);
    }

    double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
    cout << "Approximate Stock Price at T : " << avgST << endl;

    double avgVT = accumulate( Vj.begin(), Vj.end(), 0.0)/Vj.size();
    cout << "Euler Asset-or-nothing Call Price: " << exp(-optionInfo->getR()*(optionInfo->getT()-t))*avgVT << endl;
}

//void runMonteCarloSimulationForAssetOrNothing()
//{
//    unique_ptr<MonteCarlo> mc;
//    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
//    int size=100;
//    double dt=1.0/365.0;
//    double t=0.0;
//    vector<double> Sj(size,0.0),Vj(size,0.0);
//    Sj= mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(), dt,size);
//    for (int j = 0; j < Sj.size(); ++j)
//    {
//        //        cout<<Sj[j]<<endl;
//        Vj[j]=optionInfo->payOffForAssetOrNothing(Sj[j]);
//        //      cout<<Vj[j]<<endl;
//    }
//
//    double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//    cout << "Approximate Stock Price at T : " << avgST << endl;
//
//    double avgVT = accumulate( Vj.begin(), Vj.end(), 0.0)/Vj.size();
//    cout << "MC Asset-or-nothing Call Price: " << exp(-optionInfo->getR()*(optionInfo->getT()-t))*avgVT << endl;
//}

void runMonteCarloSimulationForAssetOrNothing (shared_ptr<OptionInfo> optionInfo, int MCSIM)
{
    srand(2019);
    unique_ptr<MonteCarlo> mc;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int NSIM = 5;
//    const int MCSIM=1000;
    double simSize[NSIM] = {100,200,300,400,500};
    vector<double> x(NSIM,0.0),y(NSIM,0.0),ci(NSIM,0.0);
    double dt = 1.0 / 365.0;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values = bs->calcAssetOrNothingExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(),
                                                optionInfo->getSigma(), optionInfo->getT());
//    cout<<"BS: "<<values[0]<<endl;
    for (int i = 0; i < NSIM; ++i)
    {
        double avgMCPrice=0.0;
        clock_t begin = clock();
        double MCPrice,avgSigmaHat;
        vector<double> mcPrice(MCSIM,0.0),mcError(MCSIM,0.0);
        for (int l = 0; l < MCSIM; ++l)
        {
            vector<double> Sj(simSize[i], 0.0);
            double payOff = 0.0;
            Sj = mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(), dt,
                                    simSize[i]);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
            for (int j = 0; j < Sj.size(); ++j)
            {
                payOff += optionInfo->payOffForAssetOrNothing(Sj[j]);

            }
            double avgVT = payOff / Sj.size();
//            cout<<"avgVT: "<<avgVT<<endl;
//            double sigmaHatSqr = 0.0;
//            for (int k = 0; k < Sj.size(); ++k)
//            {
//                sigmaHatSqr += (avgVT - optionInfo->payOff(Sj[k])) * (avgVT - optionInfo->payOff(Sj[k]));
//            }
//            avgSigmaHat = sigmaHatSqr / Sj.size();
            MCPrice=D*avgVT;
            mcPrice[l]=MCPrice;
            mcError[l]=(values[0]-MCPrice);
//            cout << MCPrice<<endl;
//            avgMCPrice += MCPrice;
            y[i]+=(values[0]-MCPrice);
//            ci[i]+=D * 1.96 * sqrt(avgSigmaHat / Sj.size());
        }
        mc->writeToFile("MCAONPRICE",mcPrice,i,"MCP");
        mc->writeToFile("MCAONERR",mcError,i,"ERR");
//        cout<<"y before: "<<y[i]<<endl;
        x[i]=sqrt(1.0/simSize[i]);
        y[i]=y[i]/MCSIM;
        avgMCPrice = accumulate(mcPrice.begin(), mcPrice.end(), 0.0) / mcPrice.size();
        double tempCI=0.0;
        for (int m = 0; m < MCSIM; ++m)
        {
            tempCI += (avgMCPrice-mcPrice[m])*(avgMCPrice-mcPrice[m]);
        }
        ci[i]=D*1.96*tempCI/MCSIM;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize[i] << " is :" << elapsed_secs << endl;
        cout << "MC Asset or Nothings Call Price: " << avgMCPrice <<", simSize: "<< simSize[i] <<" with CI: (+-) " << ci[i] << endl;

    }
    mc->writeToFile("MCERROR",x,y,ci,1);
}

void runAssetOrNothing()
{
    unique_ptr<BlackScholes> bs;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    vector<double> AssetOrNothingValues= bs->calcAssetOrNothingExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                                                           optionInfo->getT());
    cout<<"Asset-or-nothing Call Price: "<<AssetOrNothingValues[0]<<endl;
    cout<<"Asset-or-nothing Delta: "<<AssetOrNothingValues[1]<<endl;
}

void runMCAONDelta (shared_ptr<OptionInfo> optionInfo, int MCSIM)
{
    srand(2019);
    unique_ptr<MonteCarlo> mc;
//    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();

    const int NSIM = 4;
//    const int MCSIM = 1000;
    int sSize=100;
    int simSize;
    vector<double> x(NSIM,0.0),y(NSIM,0.0),ci(NSIM,0.0);
    double dt = 1.0 / 365.0;
    double dS=1.0;
    double t = 0.0;
    double D = exp(-optionInfo->getR() * (optionInfo->getT() - t));
    unique_ptr<BlackScholes> bs;
    vector<double> values = bs->calcAssetOrNothingExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(),
                                                optionInfo->getSigma(), optionInfo->getT());
//    cout<<"BS: "<<values[0]<<endl;
    for (int i = 0; i < NSIM; ++i)
    {
        simSize=(i+1)*sSize;
        double avgMCPricePlus=0.0,avgMCPriceMinus=0.0;
        clock_t begin = clock();
        double MCPlus,MCMinus,avgSigmaHat;
        vector<double> mcDelta(MCSIM,0.0),mcError(MCSIM,0.0);
        for (int l = 0; l < MCSIM; ++l)
        {
            vector<vector<double>> Sj;
            double payOffMinus = 0.0,payOffPlus=0.0;
            Sj = mc->genStockPricesForDelta(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),
                                            optionInfo->getSigma(), dt, dS,
                                            simSize);
//            double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
//            cout << "Approximate Stock Price at T : " << avgST << endl;
//            cout<<"Sj[0][0]"<<Sj[0][0]<<endl;
//            cout<<"Sj[1][0]"<<Sj[1][0]<<endl;
            for (int j = 0; j < simSize; ++j)
            {
                payOffPlus += optionInfo->payOffForAssetOrNothing(Sj[0][j]);

            }
//            cout<<"Pay Off Plus: "<<payOffPlus<<endl;
            for (int j = 0; j < simSize; ++j)
            {
                payOffMinus += optionInfo->payOffForAssetOrNothing(Sj[1][j]);

            }
//            cout<<"Pay Off Minus: "<<payOffMinus<<endl;
            double avgVTPlus = payOffPlus / simSize;
            double avgVTMinus = payOffMinus / simSize;
//            cout<<"avgVTPlus: "<<avgVTPlus<<endl;
//            cout<<"avgVTMinus: "<<avgVTMinus<<endl;
            MCPlus=D*avgVTPlus;
            MCMinus=D*avgVTMinus;
//            cout << MCPrice<<endl;
            avgMCPricePlus += MCPlus;
            avgMCPriceMinus += MCMinus;
            mcDelta[l]=(MCPlus-MCMinus)/(2*dS);
            mcError[l]=values[1]-mcDelta[l];
//            cout<<"avgPricePlus"<<avgMCPricePlus<<endl;
//            cout<<"avgPriceMinus"<<avgMCPriceMinus<<endl;
            double sigmaHatSqr = 0.0;
            for (int k = 0; k < simSize; ++k)
            {
                sigmaHatSqr += (avgVTPlus - optionInfo->payOff(Sj[0][k])) * (avgVTPlus - optionInfo->payOff(Sj[0][k]));
            }
            avgSigmaHat = sigmaHatSqr / Sj.size();
            ci[i]+=D * 1.96 * sqrt(avgSigmaHat / simSize);
        }
//        cout<<"y before: "<<y[i]<<endl;
        mc->writeToFile("MCAONDELTA",mcDelta,i,"MCP");
        mc->writeToFile("MCAONDELTAERR",mcError,i,"ERR");
        x[i]=sqrt(1.0/simSize);
        y[i]=y[i]/MCSIM;
        double avgMCDelta = accumulate(mcDelta.begin(), mcDelta.end(), 0.0) / mcDelta.size();
        double tempCI=0.0;
        for (int m = 0; m < MCSIM; ++m)
        {
            tempCI += (avgMCDelta-mcDelta[m])*(avgMCDelta-mcDelta[m]);
        }
        ci[i]=D*1.96*tempCI/MCSIM;
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for simSize: " << simSize << " is :" << elapsed_secs << endl;
        cout << "MC Delta : " << (avgMCPricePlus-avgMCPriceMinus)/(MCSIM*2*dS) <<", simSize: "<< simSize <<" with CI: (+-) " << ci[i] << endl;

    }
    mc->writeToFile("MCERROR",x,y,ci,1);
}

void runRNGUsingStatsBM()
{

    srand(1);
    int size=1000000;
    unique_ptr<double[]> bm{new double[size]()};
    unique_ptr<RNG> rng;
    int j=1;
    for (int i = 0; i < size; i=i+2)
    {
        bm[i]= rng->rngUsingStatsBM(0)[0];
        bm[i+1]= rng->rngUsingStatsBM(0)[1];
        if(i%10000==0)
        {
            cout<<"Generated: "<< i<<" random numbers. "<<endl;
//            srand(++j);
        }
    }
    rng->writeToFile("RNG",move(bm),size,1,"rng");
}

void runRNGUsingStatsBM1 (int size)
{

    srand(1);
    unique_ptr<double[]> bm{new double[size]()};
    unique_ptr<RNG> rng;
    bm=rng->rngUsingMTE(size);
    rng->writeToFile("RNG",move(bm),size,3,"rng");
}

void runRNGUsingStatsBM2()
{

    srand(1);
    int size=1000000;
    unique_ptr<double[]> bm{new double[size]()};
    unique_ptr<RNG> rng;
    bm=rng->rngUsingMTE(size);
    rng->writeToFile("RNG",move(bm),size,4,"rng");
}


void runWeakEulerStockPrice()
{

    srand(time(0));
    double dt=1.0/365.0;
    int size=2;
    unique_ptr<double[]> ST{new double[size]()};
    unique_ptr<RNG> rng;
    unique_ptr<Euler> euler;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    for (int i = 0; i < size; i++)
    {
//        Option( S0: 70, K: 65, T: 0.5, r: 0.07, sigma: 0.27
        ST[i]=euler->getStockPrice(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),dt);
    }
    rng->writeToFile("ST",move(ST),size,1,"ST");
}

void runWeakEulerMCCP()
{

//    srand(time(0));
    srand(2019);
    double dt=1.0/365.0;
    double CP=0.0, ST,CPFinal=0.0,sigmahatsqr, sigmaCalc;
    int size=10000;
    unique_ptr<double[]> cpVal{new double[size]()};
    unique_ptr<double[]> iVal{new double[size]()};
    unique_ptr<Euler> euler;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    double D=exp(-optionInfo->getR()*optionInfo->getT());
    for (int i = 1; i <= size; i++)
    {
//        Option( S0: 70, K: 65, T: 0.5, r: 0.07, sigma: 0.27
        ST=euler->getStockPrice(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),dt);
        //cout<<ST<<endl;
        CP+=optionInfo->payOff(ST);
        iVal[i-1]=optionInfo->payOff(ST);
        double tempCP=D*((double)1/i)*CP;
        cpVal[i-1]=tempCP;

        if(i%1000==0)
        {
            CPFinal=((double)1/i)*CP;
            for (int j = 0; j < i; ++j)
            {
                sigmaCalc=(iVal[j]-CPFinal)*(iVal[j]-CPFinal);
            }
            sigmahatsqr=((double)1/i)*sigmaCalc;
            double tempCP=D*((double)1/i)*CP;
            cout<<"Weak Euler Call Price : "<<tempCP <<" , Simulation Size: "<< i<<endl;
            cout<<"Confidence Interval : [ "<<D*(-(1.96)*sqrt((double)sigmahatsqr/i))<<" , "<<D*((1.96)*sqrt((double)sigmahatsqr/i))<<" ]" <<endl;
        }
    }
    CPFinal=((double)1/size)*CP;
    for (int j = 0; j < size; ++j)
    {
        sigmaCalc=(iVal[j]-CPFinal)*(iVal[j]-CPFinal);
    }
    sigmahatsqr=((double)1/size)*sigmaCalc;
    CP=D*((double)1/size)*CP;
    euler->writeToFile("CP",move(cpVal),size,1);
    cout<<"Weak Euler Call Price : "<<CP<<endl;
//    cout<<"Sigma hat squared: "<<CPFinal-1.96*sqrt(sigmahatsqr/size)<<endl;
    cout<<"Confidence Interval : [ "<<D*(-(1.96)*sqrt((double)sigmahatsqr/size))<<" , "<<D*((1.96)*sqrt((double)sigmahatsqr/size))<<" ]" <<endl;
}

void runWeakEulerConvergenceMCError()
{

//    srand(time(0));
    srand(1);
    const int csize=10;
    double dt=1.0/365.0;
    double CP=0.0, ST,CPFinal=0.0,sigmahatsqr, sigmaCalc;
    int size=10000;
    unique_ptr<double[]> x{new double[csize]()};
    unique_ptr<double[]> y{new double[csize]()};
    unique_ptr<double[]> ci{new double[csize]()};
//    unique_ptr<double[]> cpVal{new double[size]()};
    unique_ptr<double[]> iVal{new double[size]()};
    unique_ptr<Euler> euler;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    unique_ptr<BlackScholes> bs;
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());
    double D=exp(-optionInfo->getR()*optionInfo->getT());
    int j=0;
    for (int i = 1; i <= size; i++)
    {
        ST=euler->getStockPrice(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),dt);
        CP+=optionInfo->payOff(ST);
        iVal[i-1]=optionInfo->payOff(ST);
//        CP+=(ST-optionInfo->getK());
//        iVal[i-1]=(ST-optionInfo->getK());
        if(i%1000==0)
        {
            CPFinal=((double)1/i)*CP;
            for (int j = 0; j < i; ++j)
            {
                sigmaCalc+=(iVal[j]-CPFinal)*(iVal[j]-CPFinal);
            }
            sigmahatsqr=((double)1/i)*sigmaCalc;
            double tempCP=D*((double)1/i)*CP;
            double CI=D*((1.96)*sqrt((double)sigmahatsqr/i));
            cout<<"Weak Euler Call Price : "<<tempCP <<" CI: (+-) "<<CI<<" , Simulation Size: "<< i<<endl;
            x[j]=i;
            y[j]=values[0]-tempCP;
            ci[j]=CI;
            j++;
        }
    }
    euler->writeToFile("CIvsM",move(x),move(y),move(ci),csize,1);
}


void runWeakEulerConvergenceWithDt()
{

//    srand(time(0));
    srand(1);
    const int size=4;
    const int NSIM=10000;
    double dt[size]={0.25,0.01,0.0025,0.001};
    double ST,CPFinal=0.0,sigmahatsqr;
    unique_ptr<BlackScholes> bs;
    unique_ptr<Euler> euler;
    unique_ptr<double[]> x{new double[size]()};
    unique_ptr<double[]> y{new double[size]()};
    unique_ptr<double[]> ci{new double[size]()};
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());
    double D=exp(-optionInfo->getR()*optionInfo->getT());
    int j=0;
    for (int d = 0; d < size ; d++) //k = simulation size
    {
        double CP=0.0;
        unique_ptr<double[]> cpVal{new double[NSIM]()};
        unique_ptr<double[]> iVal{new double[NSIM]()};
        clock_t begin=clock();
        for (int i = 0; i < NSIM; i++)
        {
            ST = euler->getStockPrice(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),
                                      optionInfo->getSigma(), dt[d]);
            CP += optionInfo->payOff(ST);
            iVal[i] = optionInfo->payOff(ST);
            double tempCP = D * ((double) 1 / i) * CP;
            cpVal[i] = tempCP;
        }
        CPFinal = ((double) 1 / NSIM) * CP;
        double sigmaCalc=0;
        for (int j = 0; j < NSIM; ++j)
        {
            sigmaCalc += (iVal[j] - CPFinal) * (iVal[j] - CPFinal);
        }
        clock_t end=clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout<<"Time Taken NSIM: "<<elapsed_secs<<endl;
        sigmahatsqr = ((double) 1 / NSIM) * sigmaCalc;
        double eulerCallPrice = D * ((double) 1 / NSIM) * CP;
        double CI=D * ((1.96) * sqrt((double) sigmahatsqr / NSIM));
        cout << "Weak Euler Call Price : " << eulerCallPrice << " (+-) CI: "
             << CI << " , Simulation Size: " << NSIM << ", dt: " << dt[d]
             << endl;
        x[j]=dt[d];
        y[j]=values[0]-eulerCallPrice;
        ci[j]=CI;
        j++;
    }


    //    CPFinal=((double)1/size)*CP;
//    for (int j = 0; j < size; ++j)
//    {
//        sigmaCalc=(iVal[j]-CPFinal)*(iVal[j]-CPFinal);
//    }
//    sigmahatsqr=((double)1/size)*sigmaCalc;
//    CP=D*((double)1/size)*CP;
        euler->writeToFile("CIvsDT",move(x),move(y),move(ci),size,1);
//    cout<<"Weak Euler Call Price : "<<CP<<endl;
////    cout<<"Sigma hat squared: "<<CPFinal-1.96*sqrt(sigmahatsqr/size)<<endl;
//    cout<<"Confidence Interval : [ "<<D*(-(1.96)*sqrt((double)sigmahatsqr/size))<<" , "<<D*((1.96)*sqrt((double)sigmahatsqr/size))<<" ]" <<endl;
}


void runMCUsingStatsBM()
{

    unique_ptr<MonteCarlo> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=10000;
    double dt=1.0/365.0;
    unique_ptr<double[]> cpVal{new double[size]()};
    double CP=0.0, CPFinal=0.0, ST;
    double D=exp(-optionInfo->getR()*optionInfo->getT());
    clock_t begin=clock();
    for (int i = 1; i <= size; ++i)
    {
        ST= mc->genStockPrices3(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(), dt,
                                size);
        CP+=optionInfo->payOff(ST);
        CPFinal+=optionInfo->payOff(ST);
        double tempCP=D*((double)1/i)*CP;
        cpVal[i-1]=tempCP;
    }
    clock_t end=clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"Total Time: "<< elapsed_secs<<" seconds"<<endl;
    cout<<"MC Call Price: "<<D*((double)1/size)*CPFinal<<endl;
    mc->writeToFile("MCCP",move(cpVal),size,1);

}

void runRNG (int size)
{
    srand(1);
    vector<double> v;
    unique_ptr<RNG> rng;
     v= rng->rngUsingBM(size);
    rng->writeToFile("RNG",v,2,"rng");

}

void runWeakEulerConvergenceWithDtF()
{

//    srand(time(0));
    srand(1);

}