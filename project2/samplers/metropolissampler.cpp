#include "metropolissampler.h"
#include <ctime>
#include <random>
#include <iostream>

using std::cout;
using std::endl;

MetropolisSampler::MetropolisSampler()
{

}

MetropolisSampler::~MetropolisSampler()
{

}

double MetropolisSampler::Ratio(double ** rPosNew, double ** rPosOld, int i, double newWF, double oldWF)
{
    cout << "Overhead class default not implemented. Exiting..." << endl;
    exit(1);
    return 0.0;
}

bool MetropolisSampler::move(double ** rPosNew, double ** rPosOld, int i, double newWF, double oldWF)
{
    cout << "Overhead class default not implemented. Exiting..." << endl;
    exit(1);
    return false;
}

double MetropolisSampler::nextStep(double ** rPosOld, int i, int j)
{
    cout << "Overhead class default not implemented. Exiting..." << endl;
    exit(1);
    return 0.0;
}

double MetropolisSampler::initializePosition()
{
    cout << "Overhead class default not implemented. Exiting..." << endl;
    return 0.0;
}
