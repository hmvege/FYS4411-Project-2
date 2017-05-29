#ifndef VMC_H
#define VMC_H

#include <random>
#include "wavefunctions/wavefunctions.h"
#include "samplers/metropolissampler.h"

class VMC
{
private:
    int nParticles;
    int nDimensions;
    int acceptanceCounter = 0;
    unsigned int MCCycles;

    // Variables used globally by program
    double oldWF;           // Old wavefunction
    double newWF;           // New wavefunction
    double **rOld;          // Old positions
    double **rNew;          // New positions
    double E = 0;           // Energy
    double ESum = 0;        // Energy sum
    double ESumSquared = 0; // Energy squared sum

    void resetVariables();
    void updateParticle(int i);
    void runMetropolisStep();
    void runSDStep();
    void sampleSystem();
    void sampleSystemSD();
    void statistics(int cycles);

    WaveFunctions *WF = nullptr;
    MetropolisSampler *R = nullptr;
public:
    VMC(int new_nParticles, int new_nDimensions);
    ~VMC();
    void runVMC(unsigned int newMCCycles, unsigned int optimizationCycles, int maxSteepestDescentIterations);
    void printResults();

    // TEMP ==========================================================================================
    void diagnostics();
//    MetropolisSampler *SDR = nullptr;
    // ===============================================================================================

    // Setters
    void setWaveFunction(WaveFunctions *newWF) { WF = newWF; }
    void setNParticles(int new_nParticles) { nParticles = new_nParticles; }
    void setNDimensions(int new_nDimensions) { nDimensions= new_nDimensions; }
    void setMetropolisSampler(MetropolisSampler *newRatio) { R = newRatio; }
};

#endif // VMC_H
