#include <iostream>
#include <memory>
#include "OptionInfo.hpp"
#include "RNG.hpp"
#include "BlackScholes.hpp"
#include "EulerMethod.hpp"
#include <iomanip>



void initOptionInfo()
{
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    cout<<  optionInfo << endl;
}


void runEuler()
{
    unique_ptr<Euler> mc;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    int size=50000;
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


void runBlackScholes()
{
    unique_ptr<BlackScholes> bs;
    shared_ptr<OptionInfo> optionInfo=make_shared<OptionInfo>();
    vector<double> values= bs->calcExactValues(optionInfo->getS0(), optionInfo->getK(), optionInfo->getR(), optionInfo->getSigma(),
                                               optionInfo->getT());
    cout<<"BS Call Price: "<<values[0]<<endl;


}

int main (){

    initOptionInfo();
    runEuler();
    runBlackScholes();
    return 0;

}