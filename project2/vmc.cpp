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

void VMC::updateParticle(double **rOld, double **rNew, double &oldWF, double &newWF, int i)
{
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
    for (int i = 0; i < nParticles; i++)
    {
        updateParticle(rOld, rNew, oldWF, newWF, i);
        sampleSystem(rOld, oldWF, newWF);
    }
}

void VMC::runVMC(unsigned int newMCCycles, unsigned int optimizationCycles)
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
    R->initializePositions(rPositionsOld, rPositionsNew);
    oldWaveFunction = WF->calculate(rPositionsOld);
    while (SDCounter < maxSDIterations)
    {
        for (unsigned int i = 0; i < optimizationCycles; i++)
        {
            runMetropolisStep(rPositionsOld,rPositionsNew,oldWaveFunction,newWaveFunction);
        }
        WF->steepestDescent(rPositionsOld);
        WF->printVariationalParameters();
        SDCounter++;
    }
    // Main part of Metropolis ====================================================================
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

void VMC::sampleSystem(double ** r, double newWF, double oldWF)
{
    /*
     * Statistics that should be sampled for each run. CHANGE TO STORE IN ARRAY?
     */
    E = WF->localEnergy(r);
    ESum += E;
    ESumSquared += E*E;

//    // For steepest descent
//    dPsiAlpha       = -0.5*omega*(r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]);
//    double r12      = sqrt((r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1])); // sqrt((x1-x2)^2 + (y1-y2)^2)
//    dPsiBeta        = - a/( (beta + 1/r12)*(beta + 1/r12) );

//    dPsiBetaSum     += - a/( (beta + 1/r12)*(beta + 1/r12) );
//    dPsiAlphaSum    += -0.5*omega*(r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]);
//    dPsiEAlphaSum   += dPsiAlpha*E;
//    dPsiEBetaSum    += dPsiBeta*E;
}

void VMC::getStatistics()
{
    /*
     * Gets basic statistics of the calculations
     */
    ESum /= double(nParticles*MCCycles);        // Getting energy per particle
    ESumSquared /= double(nParticles*MCCycles);
    cout << "Energy = " << ESum << endl;
    cout << "Variance = " << (ESumSquared - ESum*ESum)/double(MCCycles) << endl;
    cout << "Acceptance rate = " << double(acceptanceCounter) / double(nParticles * MCCycles)<< endl;
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
