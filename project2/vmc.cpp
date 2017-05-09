#include <random>
#include <iostream>
#include <ctime>
#include "vmc.h"
#include "samplers/metropolissampler.h"

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
    // Setting MC cycles
    MCCycles = newMCCycles;
    // Initializing variables and configurations
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
//            rPositionsOld[i][j] = stepLength * uniform_distribution(generator) - 0.5; // OLD METHOD - WHY CENTERING AROUND 0?
            rPositionsOld[i][j] = 0.0;
//            rPositionsOld[i][j] = R->nextStep(rPositionsOld,i, j);
//            rPositionsNew[i][j] = rPositionsOld[i][j];
        }
    }
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j ++)
        {
            rPositionsOld[i][j] = R->initializePosition();
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
//                rPositionsNew[i][j] = rPositionsOld[i][j] + stepLength * uniform_distribution(generator) - 0.5; // OLD; ADD SPECIFIC SAMPLING METHOD CHOICE
//                rPositionsNew[i][j] = rPositionsOld[i][j] + stepLength * uniform_distribution(generator); // ADD SPECIFIC SAMPLING METHOD CHOICE
                rPositionsNew[i][j] = rPositionsOld[i][j] + R->nextStep(rPositionsOld,i,j);
            }
            newWaveFunction = WF->calculate(rPositionsNew); // Find the position with updated wavefunctions
//            if (uniform_distribution(generator) <= (newWaveFunction*newWaveFunction)/(oldWaveFunction*oldWaveFunction)) // OLDOLD
//            if (uniform_distribution(generator) <= R->Ratio(rPositionsOld, rPositionsNew, newWaveFunction, oldWaveFunction)) //OLD
            if (R->move(rPositionsNew, rPositionsOld, i, newWaveFunction, oldWaveFunction))
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
