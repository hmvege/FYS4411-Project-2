#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include "vmc.h"
#include "samplers/metropolissampler.h"

using std::cout;
using std::endl;

VMC::VMC(int new_nParticles, int new_nDimensions, std::string newFilename)
{
    /*
     * Constructor for the Variational Metropolis class. Arguments:
     * new_nParticles   : number of particles
     * new_nDimensions  : number of dimensions
     */
    filename = newFilename;
    setNParticles(new_nParticles);
    setNDimensions(new_nDimensions);
    oldWF = 0;
    newWF = 0;
    rOld = new double * [nParticles];
    rNew = new double * [nParticles];
}

VMC::~VMC()
{
    for (int i = 0; i < nParticles; i++)
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
     * Arguments:
     *  i   : particle being updated
     */
    R->updatePositions(rOld, rNew, i);
    newWF = WF->calculate(rNew, i);
    if (R->move(rOld, rNew, i, newWF, oldWF))
    {
        WF->updateWF();
        for (int j = 0; j < nDimensions; j++)
        {
            rOld[i][j] = rNew[i][j];
            oldWF = newWF;
        }
        acceptanceCounter++;
    }
    else
    {
        WF->revert();
        for (int j = 0; j < nDimensions; j++)
        {
            rNew[i][j] = rOld[i][j];
        }
    }
}

void VMC::runMetropolisStep(int cycle)
{
    /*
     * Metropolis update
     */
    for (int i = 0; i < nParticles; i++)
    {
        updateParticle(i);
    }
    sampleSystem(cycle);
}

void VMC::runSDStep()
{
    /*
     * Metropolis update with steepest descent
     */
    for (int i = 0; i < nParticles; i++)
    {
        updateParticle(i);
    }
    sampleSystemSD();
}

void VMC::runVMC(unsigned int newMCCycles, unsigned int optimizationCycles, int maxSteepestDescentIterations, int newMCSamplingFrequency)
{
    /*
     * Function for running the Variational Monte Carlo calculation.
     */
    // Initializing variables and configurations ==================================================
    MCCycles = newMCCycles;
    MCSamplingFrequency = newMCSamplingFrequency;
    EArr = new double[MCSamplingFrequency];
    int SDCounter = 0; // Counter for the Steepest Descent algorithm
    // Finding the optimal values for alpha and beta ==============================================
    double EOld = 0; // For checking convergence
    while (SDCounter < maxSteepestDescentIterations) // add SD convergence criteria
    {
        resetVariables();
        R->initializePositions(rOld, rNew);
        oldWF = WF->initializeWaveFunction(rOld);
        for (unsigned int i = 0; i < optimizationCycles; i++)
        {
            runSDStep();
        }
        statistics(optimizationCycles);
        WF->steepestDescent(ESum, optimizationCycles);
        WF->printVariationalParameters(SDCounter);
        SDCounter++;
        if (std::fabs(EOld - ESum) < 1e-14){ break; } // INSERT CONVERGENCE CRITERIA FUNCTION THAT CAN ADJUST STEP-SIZE!!
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
    oldWF = WF->initializeWaveFunction(rOld);
    for (unsigned int cycle = 0; cycle < MCCycles; cycle++)
    {
        runMetropolisStep(cycle);
//        if (cycle == 100)
//        {
//            printf("Planned number of cycles reached in vmc.cpp... exiting\n");
//            exit(1);
//        }
        if (((cycle+1) % MCSamplingFrequency) == 0)
        {
            writeToFile();
        }
    }
    statistics(MCCycles);
    delete [] EArr;
}

void VMC::sampleSystem(int cycle)
{
    /*
     * Base statistics that should be sampled for each run. CHANGE TO STORE IN ARRAY?
     */
    EArr[(cycle) % MCSamplingFrequency] = E;
    E = WF->localEnergy(rOld);
    ESum += E;
    ESumSquared += E*E;
}

void VMC::writeToFile()
{
    /*
     * Writing out to file every MCSamplingFrequency.
     */
    std::ofstream file;
    file.open("output/" + filename + "_MC" + std::to_string(MCCycles) + WF->getParameterString(), std::ofstream::out | std::ofstream::binary | std::ofstream::app);
    for (int i = 0; i < MCSamplingFrequency; i++)
    {
        file.write(reinterpret_cast<char*>(&EArr[i]), sizeof(double));
    }
    file.close();
}

void VMC::sampleSystemSD()
{
    /*
     * Steepest descent system sampler
     */
    E = WF->localEnergy(rOld);
    ESum += E;
    ESumSquared += E*E;
    WF->sampleSD(rOld, E);
}

void VMC::statistics(int cycles)
{
    /*
     * Gets basic statistics of the calculations
     */
    ESum /= double(MCCycles);
    ESumSquared /= double(MCCycles);
}

void VMC::printResults()
{
    /*
     * Printing results of the VMC run.
     */
    cout << "Energy:                    " << std::setprecision(15) << ESum << endl;
    cout << "Variance:                  " << std::setprecision(15) << (ESumSquared - ESum*ESum)/double(MCCycles) << endl;
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
