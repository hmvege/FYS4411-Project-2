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
    int maxSDIterations = 100;

    // Variables used globally by progra
//    double **rOld; // Old positions
//    double **rNew; // New positions
//    double WFOld; // Old wavefunction
//    double WFNew; // New wavefunction

    // Possibly temporary, will put into arrays OR WRITE TO FILE AS BINARY?
    double E = 0;
    double ESum = 0;
    double ESumSquared = 0;

    void resetVariables();

    void updateParticle(double **rOld, double **rNew, double &oldWF, double &newWF, int i);
    void runMetropolisStep(double **rOld, double **rNew, double &oldWF, double &newWF);
    void runSDStep(double **rOld, double **rNew, double &oldWF, double &newWF);
    void sampleSystem(double **r, double newWF, double oldWF);
    void sampleSystemSD(double ** r, double newWF, double oldWF);

    WaveFunctions *WF = nullptr;
    MetropolisSampler *R = nullptr;
    void checkRatio(MetropolisSampler *newRatio);
public:
    VMC(int new_nParticles, int new_nDimensions);
    ~VMC();
    void runVMC(unsigned int newMCCycles, unsigned int optimizationCycles);
    void getStatistics();

    // TEMP
    void diagnostics(double **rOld, double **rNew, double WFOld, double WFNew);

    // Setters
    void setWaveFunction(WaveFunctions *newWF) { WF = newWF; }
    void setNParticles(int new_nParticles) { nParticles = new_nParticles; }
    void setNDimensions(int new_nDimensions) { nDimensions= new_nDimensions; }
//    void setRNGSeed(double newSeed) { seed = newSeed; }
    void setMetropolisSampler(MetropolisSampler *newRatio) { R = newRatio; }
};

#endif // VMC_H
