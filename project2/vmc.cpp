#include <random>
#include <iostream>
#include <ctime>
#include "vmc.h"
#include "ratios/MetropolisRatio.h"

using std::cout;
using std::endl;

VMC::VMC()
{

}

VMC::VMC(int new_nParticles, int new_nDimensions)
{
    setNParticles(new_nParticles);
    setNDimensions(new_nDimensions);
}

VMC::~VMC()
{

}

void VMC::sampleSystem(double ** rPositionsOld, double ** rPositionsNew, double newWF, double oldWF)
{
    /*
     * Statistics that should be sampled for each run. CHANGE TO STORE IN ARRAY?
     */
    E = WF->localEnergy(rPositionsOld);
    ESum += E;
    ESumSquared += E*E;
}

void VMC::getStatistics()
{
    ESum /= double(nParticles*MCCycles);        // Getting energy per particle
    ESumSquared /= double(nParticles*MCCycles);
    cout << "Energy = " << ESum << endl;
    cout << "Variance = " << (ESumSquared - ESum*ESum) << endl;
    cout << "Acceptance rate = " << double(acceptanceCounter) / double(nParticles * MCCycles)<< endl;
}

void VMC::runVMC(unsigned int newMCCycles)
{
    // Checking if ratio has been initialized
    checkRatio(*R);
    // Setting up random generators
    if (!seed) { seed = std::time(nullptr); }
    std::mt19937_64 generator(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uniform_distribution(0,1);
    // Checking if a user defined step length has been defined
    if (!stepLength) { stepLength = 1.31; }
    // Setting MC cycles
    MCCycles = newMCCycles;
    // Initializing variables
    double oldWaveFunction = 0;
    double newWaveFunction = 0;
    double ** rPositionsOld = new double * [nParticles];
    double ** rPositionsNew = new double * [nParticles];
    for (int i = 0; i < nParticles; i++)
    {
        rPositionsOld[i] = new double[nDimensions];
        rPositionsNew[i] = new double[nDimensions];
        for (int j = 0; j < nDimensions; j ++)
        {
//            rPositionsOld[i][j] = R->nextStep()
            rPositionsOld[i][j] = stepLength * uniform_distribution(generator) - 0.5;
            rPositionsNew[i][j] = rPositionsOld[i][j];
        }
    }
    // Main part of Metropolis
    oldWaveFunction = WF->calculate(rPositionsOld);
    for (unsigned int cycle = 0; cycle < MCCycles; cycle++)
    {
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                rPositionsNew[i][j] = rPositionsOld[i][j] + stepLength * uniform_distribution(generator) - 0.5; // ADD SPECIFIC SAMPLING METHOD CHOICE
            }
            newWaveFunction = WF->calculate(rPositionsNew); // Find the position with updated wavefunctions
//            if (uniform_distribution(generator) <= (newWaveFunction*newWaveFunction)/(oldWaveFunction*oldWaveFunction))
            if (uniform_distribution(generator) <= R->Ratio(rPositionsOld, rPositionsNew, newWaveFunction, oldWaveFunction))
            { // TURN ENTIRE IF-TEST STATEMENT INTO A BOOLEAN AND TEST INSIDE RATIOS CLASS?
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
            sampleSystem(rPositionsOld, rPositionsNew, oldWaveFunction, newWaveFunction);
        }
    }
    // De-allocating memory
    for (int i = 0; i < nDimensions; i++)
    {
        delete [] rPositionsNew[i];
        delete [] rPositionsOld[i];
    }
    delete [] rPositionsNew;
    delete [] rPositionsOld;
}

void VMC::checkRatio(MetropolisRatio newRatio) {
    /*
     * Initializes a default Metropolis ratio if none is given
     */
    if (!&newRatio)
    {
        MetropolisRatio R;
        setMetropolisRatio(&R);
    }
    else
    {
        setMetropolisRatio(&newRatio);
    }
}

