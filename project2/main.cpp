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

void run2Electron(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions, double omega, double alpha, double D, double deltat, double stepLength,
                  double seed, bool impSampling, bool coulomb, string filename, string outputFolder, int MCSamplingFrequency, int numprocs, int processRank);
void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed,
                      double SDStepLength, bool runImpSampling, bool coulomb, string filename, string outputFolder, int MCSamplingFrequency, int numprocs, int processRank);
void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed, double SDStepLength,
                   bool impSampling, bool coulomb, bool jastrow, string filename, string outputFolder, int MCSamplingFrequency, int numprocs, int processRank);

/*
 * TODO GENERAL:
 * [x] Switch to new r_if(r1,r2) function in all classes
 * [x] Clean up unused variables(e.g. n_VarParams) <- CHECK ALL WARNINGS
 * [x] Add data sampling to file for every MCCount % 10000
 * [x] Parallelize code
 * [x] Clean up main to a simpler interface for the heavy data gathering run
 * [x] Setup blocking
 * [ ] Production runs for no importance sampling 2,6,12,20 electrons
 * [ ] Production run for 8e8 parallelized and non-parallelized.
 * [x] Production runs for hard-coded cases
 * [ ] Production runs for no interaction(coulombInteraction=false)
 */


int main(int numberOfArguments, char* cmdLineArguments[])
{
    // Initializing the parallelizatio process
    int numprocs, processRank;
    MPI_Init (&numberOfArguments, &cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &processRank);

    // Constants
    unsigned int MCCycles   = 1e8; // RUN THIS WITH 1 CORE!! 8e8
    unsigned int optCycles  = 1e4;
    int MCSamplingFrequency = 1e5;
    int maxSDIterations     = 250; // 0 turns it completely off, 200 is default
    int nDimensions         = 2;
    // Values for running parallel
    int nParticles[4]       = {2,6,12,20};
    double omega[5]         = {1.0, 0.5, 0.1, 0.05, 0.01};
    double alpha[4][5][2]   = {
        {{1.0,0.7}, {0.95,0.66}, {0.95,0.66}, {0.91,0.5}, {0.91,0.5}},  // N=2
        {{1.03,0.6}, {0.93,0.6}, {0.83,0.6}, {0.84,0.6}, {0.84,0.6}},   // N=6
        {{1.10,0.56}, {0.93,0.5}, {0.83,0.4}, {0.84,0.3}, {0.84,0.22}}, // N=12
        {{1.06,0.5}, {0.94,0.43}, {0.83,0.3}, {0.84,0.22}, {0.84,0.13}} // N=20
    };
    double beta[4][5]       = {
        {0.4, 0.36, 0.23, 0.2, 0.2},    // N=2
        {0.47, 0.41, 0.2, 0.15, 0.1},   // N=6
        {0.47, 0.41, 0.2, 0.15, 0.1},   // N=12
        {0.47, 0.41, 0.2, 0.15, 0.1}    // N=20
    };
    // Global setings
    std::string outputFolder= "output/no_imp";
    double D                = 0.5; // equals 0.5 in atomic units
    double deltat           = 0.001; // should be either 0.01-0.001
    double SDStepLength     = 0.001; // Steepest descent step length
//    double seed             = -10-processRank;//std::time(nullptr)-processRank;
    double seed             = std::time(nullptr)-processRank;
    bool importanceSampling = false;
    bool coulombInteraction = true;
    // Timers
    clock_t programStart, programEnd;
    clock_t runStart, runEnd;
    programStart = clock();

//    // Main loop for all different cases
//    for (int i = 0; i < 4; i++) // Default is i=0, i < 4, particles
//    {
//        if (i == 2) {
//            MCCycles = 1e7;
//        } else if (i == 3) {
//            MCCycles = 1e6;
//        }
//        for (int j = 0; j < 5; j++) // Default is j=0; j < 5, omega values
//        {
//            for (int k = 0; k < 2; k++) // Jastrow factor, jastrow off/on, default is k=0; k < 2
//            {
//                runStart = clock();
//                runNElectrons(MCCycles, optCycles, maxSDIterations, nParticles[i], nDimensions, omega[j], alpha[i][j][1-k], beta[i][j], D, deltat,seed, SDStepLength, importanceSampling, coulombInteraction, k, "NElectron", outputFolder, MCSamplingFrequency, numprocs, processRank);
//                runEnd= clock();
//                if (processRank == 0) cout << "Run complete. Time used: " << ((runEnd - runStart)/((double)CLOCKS_PER_SEC)) << endl;
//            }
//        }
//    }

    // For testing on 2 particles
//    int nParticles_2        = 2;
//    double alpha_2plain     = 0.686717; // No jastrow, 2 electrons
//    double omega_2          = 1.0;
//    double alpha_2jas       = 1.0;//0.988559; // 2 electrons
//    double beta_2jas        = 0.4;//0.398665; // 2 electrons
//    MCCycles                = 2e8;
//    outputFolder= "output/2e_plain";
//    run2Electron(MCCycles, optCycles, maxSDIterations, nParticles_2, nDimensions, omega_2, alpha_2plain, D, deltat, 1.31, seed, importanceSampling, coulombInteraction, "2ElectronPlain", outputFolder, MCSamplingFrequency, numprocs, processRank);
//    importanceSampling      = true;
//    MCCycles                = 2e8;
//    outputFolder= "output/2e_jastrow";
//    run2eImpSampling(MCCycles, optCycles, maxSDIterations, nParticles_2, nDimensions, omega_2, alpha_2jas, 1.0, beta_2jas,D, deltat, seed, SDStepLength, importanceSampling, coulombInteraction, "2ElectronJastrov", outputFolder, MCSamplingFrequency, numprocs, processRank);
//    coulombInteraction      = true;
//    importanceSampling      = true;
//    outputFolder= "output/2e_jastrowWithCoulomb";
//    run2eImpSampling(MCCycles, optCycles, maxSDIterations, nParticles_2, nDimensions, omega_2, alpha_2jas, 1.0, beta_2jas,D, deltat, seed, SDStepLength, importanceSampling, coulombInteraction, "2ElectronJastrovWithCoulomb", outputFolder, MCSamplingFrequency, numprocs, processRank);


    // Loop for unperturbed energies
    MCCycles = 1e6;
    outputFolder= "output/no_interaction";
    importanceSampling = false;
    coulombInteraction = false;
    for (int i = 0; i < 4; i++) // Default is i=0, i < 4, particles
    {
        for (int j = 0; j < 5; j++) // Default is j=0; j < 5, omega values
        {
            for (int k = 0; k < 1; k++) // Jastrow factor, jastrow off/on, default is k=0; k < 2
            {
                runStart = clock();
                runNElectrons(MCCycles, optCycles, maxSDIterations, nParticles[i], nDimensions, omega[j], alpha[i][j][1-k], beta[i][j], D, deltat,seed, SDStepLength, importanceSampling, coulombInteraction, k, "NoInteractionNElectron", outputFolder, MCSamplingFrequency, numprocs, processRank);
                runEnd= clock();
                if (processRank == 0) cout << "Run complete. Time used: " << ((runEnd - runStart)/((double)CLOCKS_PER_SEC)) << endl;
            }
        }
    }

    MPI_Finalize();

    programEnd = clock();
    for (int i = 0; i < 1e2; i++) if (processRank == 0) cout << "=";
    if (processRank == 0) cout << endl; // Printing a line
    if (processRank == 0) cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}

void printRunInfo(int MCCycles, int maxSDIterations, int nParticles, int importanceSampling,
                  int coulombInteraction, bool jastrov, double omega, double alpha, double beta, int numprocs)
{
    /*
     * Small function for printing out the configuration a run
     */
    for (int i = 0; i < 1e2; i++) { cout << "="; } cout << endl; // Printing a line
    cout << "N electrons                " << nParticles << endl;
    cout << "MC-cycles:                 " << MCCycles*numprocs << endl;
    cout << "Omega:                     " << omega << endl;
    cout << "Alpha:                     " << alpha << endl;
    if (beta != 0 && jastrov) cout << "Beta:                      " << beta << endl;
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

void run2Electron(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions, double omega,
                  double alpha, double D, double deltat, double stepLength, double seed, bool impSampling, bool coulomb,
                  std::string filename, std::string outputFolder, int MCSamplingFrequency, int numprocs, int processRank)
{
    /*
     * Function for running the two electron case.
     */
    if (processRank == 0) { printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, false, omega, alpha, 0, numprocs); }
    twoElectronPlain WF_2Electron(nParticles, nDimensions, numprocs, processRank, omega, alpha);
    WF_2Electron.setCoulombInteraction(coulomb);
    VMC VMC_2Electron(nParticles, nDimensions, outputFolder, filename, numprocs, processRank);
    VMC_2Electron.setWaveFunction(&WF_2Electron);

    if (impSampling)
    {
        ImportanceSampler importanceSampling(nParticles, nDimensions, numprocs, processRank);
        importanceSampling.initializeSampling(deltat, seed, D);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
        VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD, MCSamplingFrequency);
        VMC_2Electron.printResults();
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, numprocs, processRank);
        uniformSampling.initializeSampling(stepLength, seed);
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
        VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD, MCSamplingFrequency);
        VMC_2Electron.printResults();
    }
}

void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed, double SDStepLength,
                      bool impSampling, bool coulomb, std::string filename, std::string outputFolder, int MCSamplingFrequency, int numprocs, int processRank)
{
    /*
     * Function for running the two electron case with Jastrov factor and importance sampling.
     */
    if (processRank == 0) { printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, true, omega, alpha, beta, numprocs); }
    twoElectronJastrov WF_2Jastrov(nParticles, nDimensions, numprocs, processRank, omega, alpha, a, beta);
    WF_2Jastrov.setSDStepLength(SDStepLength);
    VMC VMC_2Electron(nParticles, nDimensions, outputFolder, filename, numprocs, processRank);
    WF_2Jastrov.setCoulombInteraction(coulomb);
    if (impSampling) {
        VMC_2Electron.setWaveFunction(&WF_2Jastrov);
        ImportanceSampler importanceSampling(nParticles, nDimensions, numprocs, processRank);
        importanceSampling.initializeSampling(deltat, seed, D);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
        VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD, MCSamplingFrequency);
        VMC_2Electron.printResults();
    }
    else
    {
        VMC_2Electron.setWaveFunction(&WF_2Jastrov);
        UniformSampling uniformSampling(nParticles, nDimensions, numprocs, processRank);
        uniformSampling.initializeSampling(1.14, seed);
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
        VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD, MCSamplingFrequency);
        VMC_2Electron.printResults();
    }
}

void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed, double SDStepLength,
                   bool impSampling, bool coulomb, bool jastrow, std::string filename, std::string outputFolder, int MCSamplingFrequency, int numprocs, int processRank)
{
    if (!jastrow && !coulomb) alpha = 1.0;
    if (processRank == 0) { printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, jastrow, omega, alpha, beta, numprocs); }
    NElectron WF_NElectron(nParticles, nDimensions, numprocs, processRank, omega, alpha, beta);
    WF_NElectron.setSDStepLength(SDStepLength);
    VMC VMC_NElectron(nParticles, nDimensions, outputFolder, filename, numprocs, processRank);
    WF_NElectron.setJastrow(jastrow);
    WF_NElectron.setCoulombInteraction(coulomb);
    VMC_NElectron.setWaveFunction(&WF_NElectron);
    if (impSampling) {
        ImportanceSampler importanceSampling(nParticles, nDimensions, numprocs, processRank);
        importanceSampling.initializeSampling(deltat, seed, D);
        VMC_NElectron.setMetropolisSampler(&importanceSampling);
        VMC_NElectron.runVMC(MCCycles,optCycles,maxNSD, MCSamplingFrequency);
        VMC_NElectron.printResults();
    }
    else
    {
        double stepLength = 1.20;
        if (omega == 1.0) {
            stepLength = 1.20;
        } else if (omega == 0.5) {
            stepLength = 1.75;
        } else if (omega == 0.1) {
            stepLength = 3.8;
        } else if (omega == 0.05) {
            stepLength = 5.2;
        } else {
            stepLength = 11.0;
        }
        UniformSampling uniformSampling(nParticles, nDimensions, numprocs, processRank);
        uniformSampling.initializeSampling(stepLength, seed);
        VMC_NElectron.setMetropolisSampler(&uniformSampling);
        VMC_NElectron.runVMC(MCCycles,optCycles,maxNSD, MCSamplingFrequency);
        VMC_NElectron.printResults();
    }
}
