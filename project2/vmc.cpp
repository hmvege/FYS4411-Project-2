#include <random>
#include <iostream>
#include <ctime>
#include "vmc.h"
#include "samplers/metropolissampler.h"

using std::cout;
using std::endl;

VMC::VMC(int new_nParticles, int new_nDimensions)
{
    /*
     * Constructor for the Variational Metropolis class. Arguments:
     * new_nParticles   : number of particles
     * new_nDimensions  : number of dimensions
     */
    setNParticles(new_nParticles);
    setNDimensions(new_nDimensions);
}

VMC::~VMC()
{

}

void VMC::updateParticle(double **rOld, double **rNew, double &oldWF, double &newWF, int i)
{
    /*
     * Function that performs a Metropolis update on a single particle
     */
    R->updatePositions(rOld, rNew, i);
    newWF = WF->calculate(rNew);
    if (R->move(rNew, rOld, i, newWF, oldWF))
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rOld[i][j] = rNew[i][j];
            oldWF = newWF;
        }
        acceptanceCounter++;
    }
    else
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew[i][j] = rOld[i][j];
        }
    }
}

void VMC::runMetropolisStep(double **rOld, double **rNew, double &oldWF, double &newWF)
{
    /*
     * Metropolis update
     */
    for (int i = 0; i < nParticles; i++)
    {
        updateParticle(rOld, rNew, oldWF, newWF, i);
        sampleSystem(rOld);
    }
}

void VMC::runSDStep(double **rOld, double **rNew, double &oldWF, double &newWF)
{
    /*
     * Metropolis update with steepest descent
     */
    for (int i = 0; i < nParticles; i++)
    {
        updateParticle(rOld, rNew, oldWF, newWF, i);
        sampleSystemSD(rOld);
    }
}

void VMC::runVMC(unsigned int newMCCycles, unsigned int optimizationCycles, int maxSteepestDescentIterations)
{
    /*
     * Function for running the Variational Monte Carlo calculation.
     */
    // Initializing variables and configurations ==================================================
    MCCycles = newMCCycles;
    int SDCounter = 0; // Counter for the Steepest Descent algorithm
    double oldWaveFunction = 0;
    double newWaveFunction = 0;
    double ** rPositionsOld = new double * [nParticles];
    double ** rPositionsNew = new double * [nParticles];
    // Finding the optimal values for alpha and beta ==============================================
    WF->printVariationalParameters();
    while (SDCounter < maxSteepestDescentIterations)
    {
        resetVariables(); // Resets E and ESum
        R->initializePositions(rPositionsOld, rPositionsNew);
        oldWaveFunction = WF->calculate(rPositionsOld);
        for (unsigned int i = 0; i < optimizationCycles; i++)
        {
            runSDStep(rPositionsOld,rPositionsNew,oldWaveFunction,newWaveFunction);
        }
        WF->steepestDescent(ESum, optimizationCycles);
        WF->printVariationalParameters();
        SDCounter++;
    }
    // Main part of Metropolis ====================================================================
    acceptanceCounter = 0;
    R->initializePositions(rPositionsOld, rPositionsNew);
    oldWaveFunction = WF->calculate(rPositionsOld);
    for (unsigned int cycle = 0; cycle < MCCycles; cycle++)
    {
        runMetropolisStep(rPositionsOld,rPositionsNew,oldWaveFunction,newWaveFunction);
    }
    // De-allocating memory =======================================================================
    for (int i = 0; i < nDimensions; i++)
    {
        delete [] rPositionsNew[i];
        delete [] rPositionsOld[i];
    }
    delete [] rPositionsNew;
    delete [] rPositionsOld;
}

void VMC::sampleSystem(double ** r)
{
    /*
     * Base statistics that should be sampled for each run. CHANGE TO STORE IN ARRAY?
     */
    E = WF->localEnergy(r);
    ESum += E;
    ESumSquared += E*E;
}

void VMC::sampleSystemSD(double ** r)
{
    /*
     * Steepest descent system sampler
     */
    sampleSystem(r);
    WF->sampleSD(r, E);
}

void VMC::getStatistics()
{
    /*
     * Gets basic statistics of the calculations
     */
    ESum /= double(nParticles*MCCycles);
    ESumSquared /= double(nParticles*MCCycles);
    cout << "Energy = " << ESum << endl;
    cout << "Variance = " << (ESumSquared - ESum*ESum)/double(MCCycles) << endl;
    cout << "Acceptance rate = " << double(acceptanceCounter) / double(nParticles * MCCycles)<< endl;
}

void VMC::resetVariables()
{
    E = 0;
    ESum = 0;
}

void VMC::diagnostics(double **rOld, double **rNew, double WFOld, double WFNew) { // REMOVE WHEN COMPLETE OR CLEAN UP!!
    cout << "WFOld = " << WFOld << endl;
    cout << "Printing rOld:" << endl;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            cout << rOld[i][j] << " ";
        }
        cout << endl;
    }
    cout << "WFNew = " << WFNew << endl;
    cout << "Printing rNew:" << endl;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            cout << rNew[i][j] << " ";
        }
        cout << endl;
    }
    R->printQMForces();
}
