#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
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
    oldWF = 0;
    newWF = 0;
    rOld = new double * [nParticles];
    rNew = new double * [nParticles];

}

VMC::~VMC()
{
    for (int i = 0; i < nDimensions; i++)
    {
        delete [] rNew[i];
        delete [] rOld[i];
    }
    delete [] rNew;
    delete [] rOld;
}

void VMC::updateParticle(int i)
{
    /*
     * Function that performs a Metropolis update on a single particle
     */
    R->updatePositions(rOld, rNew, i);
    newWF = WF->calculate(rNew);
    if (R->move(rOld, rNew, i, newWF, oldWF))
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

void VMC::runMetropolisStep()
{
    /*
     * Metropolis update
     */
    for (int i = 0; i < nParticles; i++)
    {
        updateParticle(i);
        sampleSystem();
    }
}

void VMC::runSDStep()
{
    /*
     * Metropolis update with steepest descent
     */
    for (int i = 0; i < nParticles; i++)
    {
        updateParticle(i);
//        cout << "Oops: remember to change to SD step instead" << endl; exit(1);

//        // TEMP ==========================================================================================
//        SDR->updatePositions(rOld, rNew, i);
//        newWF = WF->calculate(rNew);
//        if (SDR->move(rOld, rNew, i, newWF, oldWF))
//        {
//            for (int j = 0; j < nDimensions; j++)
//            {
//                rOld[i][j] = rNew[i][j];
//                oldWF = newWF;
//            }
//            acceptanceCounter++;
//        }
//        else
//        {
//            for (int j = 0; j < nDimensions; j++)
//            {
//                rNew[i][j] = rOld[i][j];
//            }
//        }
//        // ===============================================================================================

        sampleSystemSD();
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
    // Finding the optimal values for alpha and beta ==============================================
//    WF->printVariationalParameters();

    double EOld = 0; // For checking convergence

    while (SDCounter < maxSteepestDescentIterations) // add SD convergence criteria
    {
        resetVariables();
        R->initializePositions(rOld, rNew);
//        SDR->initializePositions(rOld, rNew); // TEMP!
        oldWF = WF->calculate(rOld);
        for (unsigned int i = 0; i < optimizationCycles; i++)
        {
            runSDStep();
        }
        statistics(optimizationCycles);
        WF->steepestDescent(ESum, optimizationCycles);

        WF->printVariationalParameters();
//        cout << "ESum " << std::setw(15) << ESum << "  EOld " << std::setw(15) << EOld << " Diff: " << std::setw(15) << std::fabs(EOld - ESum) << endl;

        SDCounter++;
        if (std::fabs(EOld - ESum) < 1e-14){ break; } // Correct convergence criteria?
        EOld = ESum;
    }

    if (maxSteepestDescentIterations != 0) // Only activates if steepest descent is used
    {
        if (SDCounter==maxSteepestDescentIterations)
        {
            cout << "Warning! No convergence found after " << SDCounter << " iterations" << endl;
        }
        else
        {
            cout << "Convergence after " << SDCounter << " steepest descent iterations" << endl;
        }
    }

    // Main part of Metropolis ====================================================================
    resetVariables();
    R->initializePositions(rOld, rNew);
    oldWF = WF->calculate(rOld);
    for (unsigned int cycle = 0; cycle < MCCycles; cycle++)
    {
        runMetropolisStep();
    }
    statistics(MCCycles);
}

void VMC::sampleSystem()
{
    /*
     * Base statistics that should be sampled for each run. CHANGE TO STORE IN ARRAY?
     */
    E = WF->localEnergy(rOld);
    ESum += E;
    ESumSquared += E*E;
}

void VMC::sampleSystemSD()
{
    /*
     * Steepest descent system sampler
     */
    sampleSystem();
    WF->sampleSD(rOld, E);
}

void VMC::statistics(int cycles)
{
    /*
     * Gets basic statistics of the calculations
     */
    ESum /= double(nParticles*cycles);
    ESumSquared /= double(nParticles*cycles);
}

void VMC::printResults()
{
    /*
     * Printing results of the VMC run.
     */
    cout << "Energy:                    " << ESum << endl;
    cout << "Variance:                  " << (ESumSquared - ESum*ESum)/double(MCCycles) << endl;
    cout << "Acceptance rate:           " << double(acceptanceCounter) / double(nParticles * MCCycles) * 100 << " %" << endl;
}

void VMC::resetVariables()
{
    /*
     * Resets variables used between the steepest desent part and the VMC main run.
     */
    E = 0;
    ESum = 0;
    acceptanceCounter = 0;
}

void VMC::diagnostics()
{
    // TEMP FUNCITON
    cout << "WFOld = " << oldWF << endl;
    cout << "Printing rOld:" << endl;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            cout << rOld[i][j] << " ";
        }
        cout << endl;
    }
    cout << "WFNew = " << newWF << endl;
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
