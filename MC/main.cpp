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
//    runMonteCarloSimulation();
//    runEuler();
    runBlackScholes();
//    cout<<endl;
//    runMonteCarloSimulationForAssetOrNothing();
//    runEulerAssetOrNothing();
//    runAssetOrNothing();
    runRNGUsingStatsBM();
//    runWeakEulerStockPrice();
//    runWeakEulerMCCP();
//    runMCUsingStatsBM();
//    runRNG();
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
            cout<<"Confidence Interval : [ "<<D*(CPFinal-(1.96)*sqrt((double)sigmahatsqr/i))<<" , "<<D*(CPFinal+(1.96)*sqrt((double)sigmahatsqr/i))<<" ]" <<endl;
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
    cout<<"Confidence Interval : [ "<<D*(CPFinal-(1.96)*sqrt((double)sigmahatsqr/size))<<" , "<<D*(CPFinal+(1.96)*sqrt((double)sigmahatsqr/size))<<" ]" <<endl;
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
        ST= mc->genStockPrices(optionInfo->getS0(), optionInfo->getT(), optionInfo->getR(), optionInfo->getSigma(), dt,
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

void runRNG()
{
    int size=10000;
    unique_ptr<double[]> bm{new double[size]()};
    vector<double> v;
    unique_ptr<RNG> rng;
    for (int i = 0; i < size; i=i+2)
    {
        v= rng->rngUsingBM(size);
    }
    for (int j = 0; j < v.size(); ++j)
    {
        bm[j]=v[j];
    }
    rng->writeToFile("RNG",move(bm),size,2,"rng");

}