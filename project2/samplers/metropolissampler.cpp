#include "metropolissampler.h"
#include <ctime>
#include <random>
#include <iostream>

using std::cout;
using std::endl;

MetropolisSampler::MetropolisSampler(int new_nParticles, int new_nDimensions, WaveFunctions *newWF)
{
    /*
     * Base constructor of any of the sampling classes.
     * Arguments:
     * new_nParticles   : the total number of particles
     * new_nDimensions  : number of dimensions, should be 2.
     * newWF            : our wavefunction, an object of the WaveFunctions class
     */
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

void MetropolisSampler::printQMForces()
{
    cout << "Wrong!!" << endl;
    exit(1);
}
