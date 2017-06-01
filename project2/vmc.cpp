#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <mpi.h>
#include "vmc.h"
#include "samplers/metropolissampler.h"

using std::cout;
using std::endl;

VMC::VMC(int new_nParticles, int new_nDimensions, std::string newFilename, int new_numprocs, int new_processRank)
{
    /*
     * Constructor for the Variational Metropolis class. Arguments:
     * new_nParticles   : number of particles
     * new_nDimensions  : number of dimensions
     */
    numprocs = new_numprocs;
    processRank = new_processRank;
    filename = newFilename + std::to_string(processRank);
    setNParticles(new_nParticles);
    setNDimensions(new_nDimensions);
    oldWF = 0;
    newWF = 0;
    rOld = new double * [nParticles];
    rNew = new double * [nParticles];
    for (int i = 0; i < nParticles; i++)
    {
        rOld[i] = new double[nDimensions];
        rNew[i] = new double[nDimensions];
    }
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
//    for (int k = 0; k < nParticles; k++) {
//        if ( k != i) {
//            for (int j=0; j < nDimensions; j++) {
//                rNew[k][j] = rOld[k][j];
//            }
//        }
//    }

    R->updatePositions(rOld, rNew, i);
    newWF = WF->calculate(rNew, i);
    if (R->move(rOld, rNew, i, newWF, oldWF))
    {
        WF->updateWF();
        for (int j = 0; j < nDimensions; j++)
        {
            rOld[i][j] = rNew[i][j];
        }
        oldWF = newWF;
        acceptanceCounter++;
    }
    else
    {
        newWF = oldWF;
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
//    double EOld = 0; // For checking convergence
    while (SDCounter < maxSteepestDescentIterations) // add SD convergence criteria
    {
        R->initializePositions(rOld, rNew);
        oldWF = WF->initializeWaveFunction(rOld);
        for (unsigned int i = 0; i < optimizationCycles; i++)
        {
//            if ((processRank == 0) && (i % 1000 == 0)) {
//                cout << "i = " << std::setw(8) << i << "  oldWF = " << std::setw(12) << oldWF << "  newWF = " << std::setw(12) << newWF << endl;
//            }
            runSDStep();
        }
        statisticsSD(optimizationCycles);
        WF->steepestDescent(ESum, optimizationCycles);
//        if (processRank == 0) WF->printVariationalParameters(SDCounter);
        SDCounter++;
//        if (std::fabs(EOld - ESum) < 1e-9)
//        {
//            cout << "ok" << endl;
//            break;
//        } // INSERT CONVERGENCE CRITERIA FUNCTION THAT CAN ADJUST STEP-SIZE!!
//        EOld = ESum;
        resetVariables();
    }
    if (maxSteepestDescentIterations != 0) { WF->finalizeSD(); } // Takes the average of the WF parmaters
    if (processRank == 0) WF->printVariationalParameters(SDCounter);
    if (maxSteepestDescentIterations != 0) // Only activates if steepest descent is used
    {
        if (SDCounter==maxSteepestDescentIterations)
        {
            if (processRank == 0) cout << "Warning! No convergence found after " << SDCounter << " iterations" << endl;
        }
        else
        {
            if (processRank == 0) cout << "Convergence after " << SDCounter << " steepest descent iterations" << endl;
        }
    }
    // Main part of Metropolis ====================================================================
    R->initializePositions(rOld, rNew);
    oldWF = WF->initializeWaveFunction(rOld);
    for (unsigned int cycle = 0; cycle < MCCycles; cycle++)
    {
//        if ((processRank == 0) && (cycle % 10000 == 0)) cout << cycle << endl;
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
    WF->localEnergy(rOld, E, EKinetic, EPotential);
    EKineticSum += EKinetic;
    EKineticSquaredSum += EKinetic*EKinetic;
    EPotentialSum += EPotential;
    EPotentialSquaredSum += EPotential*EPotential;
    ESum += E;
    ESumSquared += E*E;
    EArr[(cycle) % MCSamplingFrequency] = E;
}

void VMC::writeToFile()
{
    /*
     * Writing out to file every MCSamplingFrequency.
     */
    std::ofstream file("output/" + filename + "_Particle" + std::to_string(nParticles) + "_MC" + std::to_string(MCCycles) + WF->getParameterString(), std::ofstream::binary | std::ofstream::app);
    file.write(reinterpret_cast<const char*> (EArr), MCSamplingFrequency*sizeof(double));
    file.close();
}

void VMC::sampleSystemSD()
{
    /*
     * Steepest descent system sampler
     */
    WF->localEnergy(rOld, E, EKinetic, EPotential);
    EKineticSum += EKinetic;
    EKineticSquaredSum += EKinetic*EKinetic;
    EPotentialSum += EPotential;
    EPotentialSquaredSum += EPotential*EPotential;
    ESum += E;
    ESumSquared += E*E;
    WF->sampleSD(rOld, E);
}

void VMC::statisticsSD(int cycles)
{
    /*
     * Gets basic statistics of the calculations
     */
    ESum                    /= double(cycles);
    ESumSquared             /= double(cycles);
    EKineticSum             /= double(cycles);
    EKineticSquaredSum      /= double(cycles);
    EPotentialSum           /= double(cycles);
    EPotentialSquaredSum    /= double(cycles);
}

void VMC::statistics(int cycles)
{
    /*
     * Gets basic statistics of the calculations
     */
    double ESumTemp = 0;
    double ESumSquaredTemp = 0;
    ESum /= double(cycles);
    ESumSquared /= double(cycles);
    MPI_Reduce(&ESum, &ESumTemp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ESumSquared, &ESumSquaredTemp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    ESum = ESumTemp/double(numprocs);
    ESumSquared = ESumSquaredTemp/double(numprocs);

    double EKineticSumTemp = 0;
    double EKineticSumSquaredTemp = 0;
    EKineticSum /= double(cycles);
    EKineticSquaredSum /= double(cycles);
    MPI_Reduce(&EKineticSum, &EKineticSumTemp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&EKineticSquaredSum, &EKineticSumSquaredTemp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    EKineticSum = EKineticSumTemp/double(numprocs);
    EKineticSquaredSum = EKineticSumSquaredTemp/double(numprocs);

    double EPotentialSumTemp = 0;
    double EPotentialSumSquaredTemp = 0;
    EPotentialSum /= double(cycles);
    EPotentialSquaredSum /= double(cycles);
    MPI_Reduce(&EPotentialSum, &EPotentialSumTemp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&EPotentialSquaredSum, &EPotentialSumSquaredTemp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    EPotentialSum = EPotentialSumTemp/double(numprocs);
    EPotentialSquaredSum = EPotentialSumSquaredTemp/double(numprocs);
}

void VMC::printResults()
{
    /*
     * Printing results of the VMC run.
     */
    if (processRank == 0)
    {
        cout << "Kinetic energy:            " << std::setprecision(15) << EKineticSum << endl;
        cout << "Kinetic energy variance:   " << std::setprecision(15) << (EKineticSquaredSum - EKineticSum*EKineticSum)/double(MCCycles) << endl;
        cout << "Potential energy:          " << std::setprecision(15) << EPotentialSum << endl;
        cout << "Potential energy variance: " << std::setprecision(15) << (EPotentialSquaredSum - EPotentialSum*EPotentialSum)/double(MCCycles) << endl;
        cout << "Energy:                    " << std::setprecision(15) << ESum << endl;
        cout << "Energy variance:           " << std::setprecision(15) << (ESumSquared - ESum*ESum)/double(MCCycles) << endl;
        cout << "Acceptance rate:           " << double(acceptanceCounter) / double(nParticles * MCCycles) * 100 << " %" << endl;
    }
}

void VMC::resetVariables()
{
    /*
     * Resets variables used between the steepest descent part and the VMC main run.
     */
//    for (int i = 0; i < nParticles; i++)
//    {
//        for (int j = 0; j < nDimensions; j++)
//        {
//            rOld[i][j] = 0;
//            rNew[i][j] = 0;
//        }
//    }
    newWF = 0;
    oldWF = 0;
    E = 0;
    ESum = 0;
    ESumSquared = 0;
    EKinetic = 0;
    EKineticSum = 0;
    EKineticSquaredSum = 0;
    EPotential = 0;
    EPotentialSum = 0;
    EPotentialSquaredSum = 0;
    acceptanceCounter = 0;
}

void VMC::setMetropolisSampler(MetropolisSampler *newRatio)
{
    R = newRatio;
    R->setWaveFunction(WF);
}
