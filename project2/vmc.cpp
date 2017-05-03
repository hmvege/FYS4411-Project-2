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

void VMC::update(double ** rPositionsOld,
                 double ** rPositionsNew,
                 double &oldWaveFunction,
                 double &newWaveFunction,
                 std::mt19937_64 generator,
                 std::uniform_real_distribution<double> uniform_distribution)
{
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
            acceptanceCounter++;
        }
        else
        {
            for (int j = 0; j < nDimensions; j++)
            {
                rPositionsNew[i][j] = rPositionsOld[i][j];
            }
        }
        // Sample energy
        sampleSystem(rPositionsOld);
//        E = WF->localEnergy(rPositionsOld);
//        ESum += E;
//        ESumSquared += E*E;
    }
}

void VMC::sampleSystem(double ** rPositionsOld)
{
    // Should sample stats into an array for each run
    E = WF->localEnergy(rPositionsOld);
    ESum += E;
    ESumSquared += E*E;
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

//void VMC::initializePositions()
//{

//}

void VMC::runVMC(unsigned int MCCycles)
{
    // Setting up random generators
    if (!seed) { seed = std::time(nullptr); }
    std::mt19937_64 generator(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uniform_distribution(0,1);

    double stepLength = 1.31; // Should generalize this, 1.31 gives 50% acceptance rate

    double oldWaveFunction = 0;
    double newWaveFunction = 0;

//    double E = 0;
//    double ESum = 0;
//    double ESumSquared = 0;
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
//        update(rPositionsOld, rPositionsNew, oldWaveFunction, newWaveFunction, generator, uniform_distribution);
//        sampleSystem(rPositionsOld);

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
                acceptanceCounter++;
            }
            else
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    rPositionsNew[i][j] = rPositionsOld[i][j];
                }
            }

//            // Sample energy
            sampleSystem(rPositionsOld);
//            E = WF->localEnergy(rPositionsOld);
//            ESum += E;
//            ESumSquared += E*E;
        }
    }
    ESum /= double(nParticles*MCCycles);        // Getting energy per particle
    ESumSquared /= double(nParticles*MCCycles);
    cout << "Energy = " << ESum << endl;
    cout << "Variance = " << (ESumSquared - ESum*ESum) << endl;
    cout << "Acceptance rate = " << double(acceptanceCounter) / double(nParticles * MCCycles)<< endl;


    for (int i = 0; i < nDimensions; i++)
    {
        delete [] rPositionsNew[i];
        delete [] rPositionsOld[i];
    }
    delete [] rPositionsNew;
    delete [] rPositionsOld;
}
