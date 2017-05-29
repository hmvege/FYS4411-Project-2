#include "metropolissampler.h"
#include <ctime>
#include <random>
#include <iostream>

using std::cout;
using std::endl;

MetropolisSampler::MetropolisSampler(int new_nParticles, int new_nDimensions, WaveFunctions *newWF)
{
    nParticles = new_nParticles;
    nDimensions = new_nDimensions;
    setWaveFunction(newWF);
}

MetropolisSampler::~MetropolisSampler()
{

}

double MetropolisSampler::Ratio(double ** rOld, double ** rNew, int i, double newWF, double oldWF)
{
    cout << "Overhead class default of Ratio not implemented in metropolissampler.cpp. Exiting..." << endl;
    exit(1);
    return 0.0;
}

bool MetropolisSampler::move(double ** rOld, double ** rNew, int i, double newWF, double oldWF)
{
    cout << "Overhead class default of Next not implemented in metropolissampler.cpp. Exiting..." << endl;
    exit(1);
    return false;
}

double MetropolisSampler::nextStep(double ** rOld, int i, int j)
{
    cout << "Overhead class default of nextStep not implemented in metropolissampler.cpp. Exiting..." << endl;
    exit(1);
    return 0.0;
}

void MetropolisSampler::updatePositions(double ** rOld, double ** rNew, int k)
{
    cout << "Overhead class default of updatePositions not implemented in metropolissampler.cpp. Exiting..." << endl;
    exit(1);
}

void MetropolisSampler::initializePositions(double **rOld, double **rNew)
{
    cout << "Overhead class default of initializePositsions not implemented in metropolissampler.cpp. Exiting..." << endl;
    exit(1);
}

// TMEP
void MetropolisSampler::printQMForces()
{
    cout << "Wrong!!" << endl;
    exit(1);
}
