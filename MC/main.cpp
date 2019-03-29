#include <iostream>
#include <memory>
#include "OptionInfo.hpp"
#include "RNG.hpp"
#include "MonteCarlo.hpp"
#include "EulerMethod.hpp"
#include "BlackScholes.hpp"
#include <iomanip>

void initOptionInfo();
void runMonteCarloSimulation();
void runRNG();
void runEuler();
void runBlackScholes();
void runAssetOrNothing();

int main ()
{
    initOptionInfo();
    runMonteCarloSimulation();
    runEuler();
    runBlackScholes();
    runAssetOrNothing();
    return 0;
}

void initOptionInfo()
{
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    cout<<  optionInfo<<endl;
}

void runEuler()
{
    unique_ptr<Euler> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=100;
    double t = 0.0;

    vector<double> Sj(size,0.0),Vj(size,0.0),SjMinusDT(size,0.0),SjPlusDT(size,0.0);
    vector<double> VjMinusDT(size,0.0),VjPlusDT(size,0.0);
    Sj = mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),size);
    for (int j = 0; j < Sj.size(); ++j)
    {
        Vj[j]=optionInfo->payOff(Sj[j]);

    }

    double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
    cout << "Approximate Stock Price at T : " << avgST << endl;

    double avgVT = accumulate( Vj.begin(), Vj.end(), 0.0)/Vj.size();
    cout << "Euler Call Price: " << exp(-optionInfo->getR()*(optionInfo->getT()-t))*avgVT << endl;



}

void runMonteCarloSimulation()
{
    unique_ptr<MonteCarlo> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=100;
    double dt=1.0/365.0;
    double t=0.0;
    vector<double> Sj(size,0.0),Vj(size,0.0),SjMinusDT(size,0.0),SjPlusDT(size,0.0);
    vector<double> VjMinusDT(size,0.0),VjPlusDT(size,0.0);
    Sj= mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),size);
    for (int j = 0; j < Sj.size(); ++j)
    {
    //        cout<<Sj[j]<<endl;
            Vj[j]=optionInfo->payOff(Sj[j]);

    }

    double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
    cout << "Approximate Stock Price at T : " << avgST << endl;

    double avgVT = accumulate( Vj.begin(), Vj.end(), 0.0)/Vj.size();
    cout << "MC Call Price: " << exp(-optionInfo->getR()*(optionInfo->getT()-t))*avgVT << endl;

    //DELTA SIMULATION
    SjMinusDT= mc->genStockPrices(optionInfo->getS0(), optionInfo->getT()-dt, optionInfo->getR(),optionInfo->getSigma(),size);
    SjPlusDT = mc->genStockPrices(optionInfo->getS0(), optionInfo->getT()+dt, optionInfo->getR(),optionInfo->getSigma(),size);
    for (int j = 0; j < Sj.size(); ++j)
    {
        //        cout<<Sj[j]<<endl;
        VjMinusDT[j]=optionInfo->payOff(SjMinusDT[j]);
        VjPlusDT[j]=optionInfo->payOff(SjPlusDT[j]);
    }

    double avgVTMinusDT = accumulate( VjMinusDT.begin(), VjMinusDT.end(), 0.0)/VjMinusDT.size();
    double avgVTPlusDT = accumulate( VjPlusDT.begin(), VjPlusDT.end(), 0.0)/VjPlusDT.size();
    cout << "MC Delta: " << ((exp(-optionInfo->getR()*(optionInfo->getT()+dt-t))*avgVTPlusDT)-(exp(-optionInfo->getR()*(optionInfo->getT()-dt-t))*avgVTMinusDT))/(2*dt) << endl;

}

void runBlackScholes()
{
    unique_ptr<BlackScholes> bs;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                  optionInfo->getT());
    cout<<"BS Call Price: "<<values[0]<<endl;
    cout<<"BS Delta: "<<values[1]<<endl;

}

void runAssetOrNothing()
{
    unique_ptr<BlackScholes> bs;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    vector<double> AssetOrNothingValues= bs->calcAssetOrNothingExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                                                           optionInfo->getT(),optionInfo->getDY());
    cout<<"Asset-or-nothing Call Price: "<<AssetOrNothingValues[0]<<endl;
    cout<<"Asset-or-nothing Delta: "<<AssetOrNothingValues[1]<<endl;
}