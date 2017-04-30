#ifndef VMC_H
#define VMC_H

#include "wavefunctions/wavefunctions.h"

class VMC
{
private:
    int nParticles;
    int nDimensions;
    int acceptanceCounter;
    double seed = false;
    double stepLength = 1.0;

    void update(double **oldPositions, double **newPositions);
    void sampleSystem();
    double R();

    WaveFunctions *WF = nullptr;
public:
    VMC();
    ~VMC();
    void runVMC(unsigned int MCCycles);

    void getStatistics();

    // Setters
    void setStepLength(double newStepLength) { stepLength = newStepLength; }
    void setWaveFunction(WaveFunctions *newWF) { WF = newWF; }
    void setNParticles(int new_nParticles) { nParticles = new_nParticles; }
    void setNDimensions(int new_nDimensions) { nDimensions= new_nDimensions; }
    void setRNGSeed(double newSeed) { seed = newSeed; }
};

#endif // VMC_H
