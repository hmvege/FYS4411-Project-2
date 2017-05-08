#ifndef VMC_H
#define VMC_H

#include <random>
#include "wavefunctions/wavefunctions.h"
#include "ratios/MetropolisRatio.h"

class VMC
{
private:
    int nParticles;
    int nDimensions;
    int acceptanceCounter = 0;
    int MCCycles;
    double seed = false;
    double stepLength = false;

    // Possibly temporary, will put into arrays?
    double E = 0;
    double ESum = 0;
    double ESumSquared = 0;

//    void update(double **rPositionsOld, double **rPositionsNew, double &oldWaveFunction, double &newWaveFunction);
    void sampleSystem(double **rPositionsOld, double **rPositionsNew, double newWF, double oldWF);
//    double R(double ** rPositionsOld, double ** rPositionsNew, double newWF, double oldWF);

    WaveFunctions *WF = nullptr;
    MetropolisRatio *R = nullptr;
    void checkRatio(MetropolisRatio newRatio);
public:
    VMC();
    VMC(int new_nParticles, int new_nDimensions);
    ~VMC();
    void runVMC(unsigned int newMCCycles);
    void getStatistics();

    // Setters
    void setStepLength(double newStepLength) { stepLength = newStepLength; }
    void setWaveFunction(WaveFunctions *newWF) { WF = newWF; }
    void setNParticles(int new_nParticles) { nParticles = new_nParticles; }
    void setNDimensions(int new_nDimensions) { nDimensions= new_nDimensions; }
    void setRNGSeed(double newSeed) { seed = newSeed; }
    void setMetropolisRatio(MetropolisRatio *newRatio) { R = newRatio; }
};

#endif // VMC_H
