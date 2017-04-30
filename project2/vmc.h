#ifndef VMC_H
#define VMC_H

#include <random>
#include "wavefunctions/wavefunctions.h"

class VMC
{
private:
    int nParticles;
    int nDimensions;
    int acceptanceCounter = 0;
    double seed = false;
    double stepLength = 1.0;

    // Possibly temporary, will put into arrays?
    double E;
    double ESum;
    double ESumSquared;

//    void update(double **rPositionsOld, double **rPositionsNew, double &oldWaveFunction, double &newWaveFunction);
    void update(double ** rPositionsOld,
                     double ** rPositionsNew,
                     double &oldWaveFunction,
                     double &newWaveFunction,
                     std::mt19937_64 generator,
                     std::uniform_real_distribution<double> uniform_distribution);
    void sampleSystem(double **rPositionsOld);
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
