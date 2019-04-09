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
void runRNG (int size);
void runEuler();
void runBlackScholes();
void runEulerAssetOrNothing();
void runMonteCarloSimulationForAssetOrNothing();
void runAssetOrNothing();
void runRNGUsingStatsBM();
void runRNGUsingStatsBM1 (int size);
void runWeakEulerStockPrice();
void runWeakEulerMCCP();
void runMCUsingStatsBM();
void runWeakEulerConvergenceWithDt();
void runWeakEulerConvergenceMCError();
void runRNGUsingStatsBM2();

int main ()
{
    initOptionInfo();
//    runMonteCarloSimulation();
    runBlackScholes();
    runEuler();
//    cout<<endl;
//    runMonteCarloSimulationForAssetOrNothing();
//    runEulerAssetOrNothing();
//    runAssetOrNothing();
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
//    runRNG(10000000);
//    end=clock();
//    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    cout<<"Time Taken vector: "<<elapsed_secs<<endl;
//    runWeakEulerConvergenceWithDt();
//    runWeakEulerConvergenceMCError();
//    runRNGUsingStatsBM2();
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
    srand(2019);
    unique_ptr<Euler> euler;
    shared_ptr<OptionInfo> optionInfo = make_shared<OptionInfo>();
    const int dtSize = 5;
    double dt[dtSize] = {0.5, 0.2, 0.1, 0.05,0.02};
    vector<double> x(dtSize),y(dtSize),ci(dtSize);
    int simSize = 1000;
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
        cout << "SJ Size:" << Sj.size() << endl;
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
        cout << "Euler Call Price: " << eulerCallPrice  << " with CI: (+-) " << D * 1.96 * sqrt(avgSigmaHat / Vj.size()) << endl;
        x[i]=dt[i];
        y[i]=values[0]-eulerCallPrice;
        ci[i]=D * 1.96 * sqrt(avgSigmaHat / Vj.size());
    }
    euler->writeToFile("EULER",x,y,ci,1);
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