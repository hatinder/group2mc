#include <iostream>
#include <memory>
#include "OptionInfo.hpp"
#include "RNG.hpp"
#include "MonteCarlo.hpp"
#include "EulerMethod.hpp"
#include "BlackScholes.hpp"
#include <iomanip>


random_device rdStats;

void initOptionInfo();
void runMonteCarloSimulation();
void runRNG();
void runEuler();
void runBlackScholes();
void runEulerAssetOrNothing();
void runMonteCarloSimulationForAssetOrNothing();
void runAssetOrNothing();
void runRNGUsingStatsBM();
void runWeakEulerStockPrice();
void runWeakEulerMCCP();
void runMCUsingStatsBM();

int main ()
{
    initOptionInfo();
    runMonteCarloSimulation();
    runEuler();
    runBlackScholes();
//    cout<<endl;
//    runMonteCarloSimulationForAssetOrNothing();
//    runEulerAssetOrNothing();
//    runAssetOrNothing();
//    runRNGUsingStatsBM();
//    runWeakEulerStockPrice();
//    runWeakEulerMCCP();
//    runMCUsingStatsBM();
    return 0;
}

void initOptionInfo()
{
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    cout<<  optionInfo<<endl;
//    shared_ptr<OptionInfo> optionInfo2=make_shared<OptionInfo>(70.0, 65.0, 0.5, 0.07, 0.27);
//    cout<<  optionInfo2<<endl;

}

void runEuler()
{
    unique_ptr<Euler> euler;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=100;
    double t = 0.0;

    vector<double> Sj(size,0.0),Vj(size,0.0),SjMinusDT(size,0.0),SjPlusDT(size,0.0);
    vector<double> VjMinusDT(size,0.0),VjPlusDT(size,0.0);
    Sj = euler->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),size);
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
    int size=10000;
    double dt=1.0/365.0;
    double t=0.0;
    double D=exp(-optionInfo->getR()*(optionInfo->getT()-t));
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
    cout << "MC Call Price: " << D*avgVT <<" , Discounting Factor:"<<D<< endl;

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

void runEulerAssetOrNothing()
{
    unique_ptr<Euler> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=100;
    double t = 0.0;

    vector<double> Sj(size,0.0),Vj(size,0.0);
    Sj = mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),size);
    for (int j = 0; j < Sj.size(); ++j)
    {
        Vj[j]=optionInfo->payOffForAssetOrNothing(Sj[j]);
    }

    double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
    cout << "Approximate Stock Price at T : " << avgST << endl;

    double avgVT = accumulate( Vj.begin(), Vj.end(), 0.0)/Vj.size();
    cout << "Euler Asset-or-nothing Call Price: " << exp(-optionInfo->getR()*(optionInfo->getT()-t))*avgVT << endl;
}

void runMonteCarloSimulationForAssetOrNothing()
{
    unique_ptr<MonteCarlo> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=100;
    double dt=1.0/365.0;
    double t=0.0;
    vector<double> Sj(size,0.0),Vj(size,0.0);
    Sj= mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),size);
    for (int j = 0; j < Sj.size(); ++j)
    {
        //        cout<<Sj[j]<<endl;
        Vj[j]=optionInfo->payOffForAssetOrNothing(Sj[j]);
        //      cout<<Vj[j]<<endl;
    }

    double avgST = accumulate( Sj.begin(), Sj.end(), 0.0)/Sj.size();
    cout << "Approximate Stock Price at T : " << avgST << endl;

    double avgVT = accumulate( Vj.begin(), Vj.end(), 0.0)/Vj.size();
    cout << "MC Asset-or-nothing Call Price: " << exp(-optionInfo->getR()*(optionInfo->getT()-t))*avgVT << endl;
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


void runRNGUsingStatsBM()
{

    srand(time(0));
    int size=10000;
    unique_ptr<double[]> bm{new double[size]()};
    unique_ptr<RNG> rng;
    for (int i = 0; i < size; i=i+2)
    {
        bm[i]= rng->rngUsingStatsBM(0)[0];
        bm[i+1]= rng->rngUsingStatsBM(0)[1];
    }
    rng->writeToFile("RNG",move(bm),size,1,"rng");
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

    srand(time(0));
    double dt=1.0/(2*365.0);
    double CP=0.0, ST;
    int size=10000;
    unique_ptr<double[]> cpVal{new double[size]()};
    unique_ptr<Euler> euler;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    double D=exp(-optionInfo->getR()*optionInfo->getT());
    for (int i = 0; i < size; i++)
    {
//        Option( S0: 70, K: 65, T: 0.5, r: 0.07, sigma: 0.27
        ST=euler->getStockPrice(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(),optionInfo->getSigma(),dt);
        //cout<<ST<<endl;
        CP+=optionInfo->payOff(ST);
        cpVal[i]=D*(ST-optionInfo->getK());
    }
    CP=D*((double)1/size)*CP;
    euler->writeToFile("CP",move(cpVal),size,1);
    cout<<"Weak Euler Call Price : "<<CP<<" , Discount Factor: "<<D<<endl;
}

void runMCUsingStatsBM()
{

    unique_ptr<MonteCarlo> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=100;
    double dt=1.0/365.0;
    unique_ptr<double[]> cpVal{new double[size]()};
    double CP=0.0, ST;
    double D=exp(-optionInfo->getR()*optionInfo->getT());
    for (int j = 1; j < size; ++j)
    {
        CP=0.0;
        for (int i = 0; i < j; ++i)
        {
            ST= mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(), dt,
                                   size);
            CP+=optionInfo->payOff(ST);

        }
        CP=D*((double)1/j)*CP;
        cpVal[j]=CP;
        cout<<"Monte Carlo Call Price : "<<CP<<" , size: "<<j<<endl;
    }

    mc->writeToFile("MCCP",move(cpVal),size,1);

}