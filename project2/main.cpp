#include <iostream>
#include <cmath>
#include <ctime>
#include <mpi.h>
#include "vmc.h"
#include "wavefunctions/twoelectronplain.h"
#include "wavefunctions/twoelectronjastrov.h"
#include "wavefunctions/n_electron_wf/nelectron.h"
#include "samplers/metropolissampler.h"
#include "samplers/uniformsampling.h"
#include "samplers/importancesampler.h"
#include "functions.h"

using namespace std;

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double stepLength,
                  double seed, bool impSampling, bool coulomb, string filename, int MCSamplingFrequency, int numprocs, int processRank);
void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed,
                      double SDStepLength, bool runImpSampling, bool coulomb, string filename, int MCSamplingFrequency, int numprocs, int processRank);
void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed, double SDStepLength,
                   bool impSampling, bool coulomb, bool jastrow, string filename, int MCSamplingFrequency, int numprocs, int processRank);

/*
 * TODO GENERAL:
 * [x] Switch to new r_if(r1,r2) function in all classes
 * [x] Clean up unused variables(e.g. n_VarParams) <- CHECK ALL WARNINGS
 * [x] Add data sampling to file for every MCCount % 10000
 * [x] Parallelize code
 * [x] Clean up main to a simpler interface for the heavy data gathering run
 * [ ] Setup blocking
 */


int main(int numberOfArguments, char* cmdLineArguments[])
{
    // Initializing the parallelizatio process
    int numprocs, processRank;
    MPI_Init (&numberOfArguments, &cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &processRank);

    // Constants
    unsigned int MCCycles   = 4e8/numprocs;
    unsigned int optCycles  = 1e5;
    int MCSamplingFrequency = 1e6;
    int maxSDIterations     = 200; // 0 turns it completely off
//    int nParticles          = 2;
    int nDimensions         = 2;
    // Values for running parallel
    int nParticles[4]       = {2,6,12,20};
    double omega[5]         = {1.0, 0.5, 0.1, 0.05, 0.1};
    double alpha[4][5] = {
        {1.0, 0.95, 0.95, 0.91, 0.91},
        {1.03, 0.93, 0.83, 0.84, 0.84},
        {1.10, 0.93, 0.83, 0.84, 0.84},
        {1.06, 0.94, 0.83, 0.84, 0.84}
    };
    double beta[4][5] = {
        {0.4, 0.36, 0.23, 0.2, 0.2},
        {0.47, 0.41, 0.2, 0.15, 0.1},
        {0.47, 0.41, 0.2, 0.15, 0.1},
        {0.47, 0.41, 0.2, 0.15, 0.1}
    };
//    double omega            = 1.0;
//    double alpha            = 1.0;//0.988559; // 2 electrons
//    double beta             = 0.4;//0.398665; // 2 electrons
//    double alpha            = 1.03741;    // 6 electrons
//    double beta             = 0.472513;   // 6 electrons
//    double alpha            = 1.10364;  // 12 electrons
//    double beta             = 0.468861; // 12 electrons
//    double alpha            = 0.686717; // No jastrow, 2 electrons
//    double alpha            = 0.599; // No jastrow, 6 electrons
//    double alpha            = 0.569619; // No jastrow, 12 electrons
    double D                = 0.5; // equals 0.5 in atomic units
    double deltat           = 0.001; // should be either 0.01-0.001
    double SDStepLength     = 0.001; // Steepest descent step length
    double seed             = -1-processRank;//std::time(nullptr);
    bool importanceSampling = false;
    bool coulombInteraction = true;
    bool jastrowFactor      = true;

    clock_t programStart, programEnd;
    clock_t runStart, runEnd;
    programStart = clock();

//    run2Electron(MCCycles, nParticles, nDimensions, omega, alpha, 1.31, seed, importanceSampling, coulombInteraction, "2ElectronPlain", MCSamplingFrequency, numprocs, processRank);
//    run2eImpSampling(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, 1.0, beta,D, deltat, seed, SDStepLength, importanceSampling, coulombInteraction, "2ElectronJastrov", MCSamplingFrequency, numprocs, processRank);
//    runNElectrons(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, beta, D, deltat,seed, SDStepLength, importanceSampling, coulombInteraction, jastrowFactor, "NElectron", MCSamplingFrequency, numprocs, processRank);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            runStart = clock();
            runNElectrons(MCCycles, optCycles, maxSDIterations, nParticles[i], nDimensions, omega[j], alpha[i][j], beta[i][j], D, deltat,seed, SDStepLength, importanceSampling, coulombInteraction, jastrowFactor, "NElectron", MCSamplingFrequency, numprocs, processRank);
            runEnd= clock();
            if (processRank == 0) cout << "Run complete. Time used: " << ((runEnd - runStart)/((double)CLOCKS_PER_SEC)) << endl;
        }
    }

    MPI_Finalize();

    programEnd = clock();
    for (int i = 0; i < 1e2; i++)
    {
        if (processRank == 0)
        {
            cout << "=";
        }
    }
    if (processRank == 0) cout << endl; // Printing a line
    if (processRank == 0) cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}

void printRunInfo(int MCCycles, int maxSDIterations, int nParticles, int importanceSampling,
                  int coulombInteraction, bool jastrov, double omega, double alpha, double beta)
{
    /*
     * Small function for printing out the configuration a run
     */
    for (int i = 0; i < 1e2; i++) { cout << "="; } cout << endl; // Printing a line
    cout << "N electrons                " << nParticles << endl;
    cout << "MC-cycles:                 " << MCCycles << endl;
    cout << "Omega:                     " << omega << endl;
    cout << "Alpha:                     " << alpha << endl;
    if (beta != 0) cout << "Beta:                      " << beta << endl;
    cout << "Jastrov                    ";
    if (jastrov) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "Steepest descent:          ";
    if (maxSDIterations > 0) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "Importance sampling:       ";
    if (importanceSampling) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "Coulomb interaction:       ";
    if (coulombInteraction) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "RESULTS:" << endl;
}

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha,
                  double stepLength, double seed, bool impSampling, bool coulomb, std::string filename,
                  int MCSamplingFrequency, int numprocs, int processRank)
{
    /*
     * Function for running the two electron case.
     */
    if (processRank == 0) { printRunInfo(MCCycles, 0, nParticles, impSampling, coulomb, false, omega, alpha, 0); }
    twoElectronPlain WF_2Electron(nParticles, nDimensions, omega, alpha);
    WF_2Electron.setCoulombInteraction(coulomb);
    VMC VMC_2Electron(nParticles, nDimensions, filename, numprocs, processRank);

    if (impSampling)
    {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_2Electron);
        importanceSampling.initializeSampling(0.001, seed, 0.5);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_2Electron);
        uniformSampling.initializeSampling(stepLength, seed);
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
    }
    VMC_2Electron.setWaveFunction(&WF_2Electron);
    VMC_2Electron.runVMC(MCCycles, 0, 0, MCSamplingFrequency);
    VMC_2Electron.printResults();
}

void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed, double SDStepLength,
                      bool impSampling, bool coulomb, std::string filename, int MCSamplingFrequency, int numprocs, int processRank)
{
    /*
     * Function for running the two electron case with Jastrov factor and importance sampling.
     */
    if (processRank == 0) { printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, true, omega, alpha, beta); }
    twoElectronJastrov WF_2Jastrov(nParticles, nDimensions, omega, alpha, a, beta);
    WF_2Jastrov.setSDStepLength(SDStepLength);
    VMC VMC_2Electron(nParticles, nDimensions, filename, numprocs, processRank);
    WF_2Jastrov.setCoulombInteraction(coulomb);
    if (impSampling) {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_2Jastrov);
        importanceSampling.initializeSampling(deltat, seed, D);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_2Jastrov);
        uniformSampling.initializeSampling(1.14, seed);
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
    }

    VMC_2Electron.setWaveFunction(&WF_2Jastrov);
    VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD, MCSamplingFrequency);
    VMC_2Electron.printResults();
}

void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed, double SDStepLength,
                   bool impSampling, bool coulomb, bool jastrow, std::string filename, int MCSamplingFrequency, int numprocs, int processRank)
{
    if (processRank == 0) { printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, jastrow, omega, alpha, beta); }
    NElectron WF_NElectron(nParticles, nDimensions, omega, alpha, beta);
    WF_NElectron.setSDStepLength(SDStepLength);
    VMC VMC_NElectron(nParticles, nDimensions, filename, numprocs, processRank);
    WF_NElectron.setCoulombInteraction(coulomb);
    WF_NElectron.setJastrow(jastrow);
    if (impSampling) {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_NElectron);
        importanceSampling.initializeSampling(deltat, seed, D);        
        VMC_NElectron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_NElectron);
        uniformSampling.initializeSampling(1.14, seed);
        VMC_NElectron.setMetropolisSampler(&uniformSampling);
    }
    VMC_NElectron.setWaveFunction(&WF_NElectron);
    VMC_NElectron.runVMC(MCCycles,optCycles,maxNSD, MCSamplingFrequency);
    VMC_NElectron.printResults();
}
