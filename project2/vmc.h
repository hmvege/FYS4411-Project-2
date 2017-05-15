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
    int MCCycles;

    // Possibly temporary, will put into arrays?
    double E = 0;
    double ESum = 0;
    double ESumSquared = 0;

//    void update(double **rPositionsOld, double **rPositionsNew, double &oldWaveFunction, double &newWaveFunction);
    void sampleSystem(double **rPositionsOld, double **rPositionsNew, double newWF, double oldWF);

    WaveFunctions *WF = nullptr;
    MetropolisSampler *R = nullptr;
    void checkRatio(MetropolisSampler *newRatio);
public:
    VMC();
    VMC(int new_nParticles, int new_nDimensions);
    ~VMC();
    void runVMC(unsigned int newMCCycles);
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
