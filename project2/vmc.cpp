#include <random>
#include <iostream>
#include <ctime>
#include "vmc.h"

using std::cout;
using std::endl;

VMC::VMC()
{

}

VMC::~VMC()
{

}

void VMC::update(double ** oldPositions, double ** newPositions) // Finish this
{
    double oldWaveFunction = 0;
    double newWaveFunction = 0;
}

void VMC::sampleSystem()
{
    // Should sample stats into an array for each run

}

double VMC::R()
{
    /*
     * VMC ratio
     */
    return 0.0;
}

void VMC::getStatistics()
{

}

void VMC::runVMC(unsigned int MCCycles)
{
    // Setting up random generators
    if (!seed) { seed = std::time(nullptr); }
    std::mt19937_64 generator(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uniform_distribution(0,1);

    double stepLength = 1.0; // Should generalize this

    // ADD ACCEPT-REJECT COUNTING

    double oldWaveFunction = 0;
    double newWaveFunction = 0;

    double E = 0;
    double ESum = 0;
    double ESumSquared = 0;
    double * MCSamples = new double[MCCycles];
    double ** rPositionsOld = new double * [nParticles];
    double ** rPositionsNew = new double * [nParticles];

    for (int i = 0; i < nParticles; i++)
    {
        rPositionsOld[i] = new double[nDimensions];
        rPositionsNew[i] = new double[nDimensions];
        for (int j = 0; j < nDimensions; j ++)
        {
            rPositionsOld[i][j] = stepLength * uniform_distribution(generator) - 0.5;
            rPositionsNew[i][j] = rPositionsOld[i][j];
        }
    }

    oldWaveFunction = WF->calculate(rPositionsOld);
    for (unsigned int cycle = 0; cycle < MCCycles; cycle++)
    {
//        MCSamples[i] = 0;
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                rPositionsNew[i][j] = rPositionsOld[i][j] + stepLength * uniform_distribution(generator) - 0.5;
            }
            newWaveFunction = WF->calculate(rPositionsNew); // Find the position with updated wavefunctions
            if (uniform_distribution(generator) <= (newWaveFunction*newWaveFunction)/(oldWaveFunction*oldWaveFunction))
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    rPositionsOld[i][j] = rPositionsNew[i][j];
                    oldWaveFunction = newWaveFunction;
                }
            }
            else
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    rPositionsNew[i][j] = rPositionsOld[i][j];
                }
            }
            E = WF->localEnergy(rPositionsOld);
            ESum += E;
            ESumSquared += E*E;
        }
    }
    ESum /= double(nParticles*MCCycles);        // Getting energy per particle
    ESumSquared /= double(nParticles*MCCycles);
    cout << "Energy = " << ESum << endl;
    cout << "Variance = " << (ESumSquared - ESum*ESum) << endl;


    for (int i = 0; i < nDimensions; i++)
    {
        delete [] rPositionsNew[i];
        delete [] rPositionsOld[i];
    }
    delete [] MCSamples; // Fix this!! Gives warning!!
    delete [] rPositionsNew;
    delete [] rPositionsOld;
}
