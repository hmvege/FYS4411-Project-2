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
    int maxSDIterations = 0;

    // Possibly temporary, will put into arrays OR WRITE TO FILE AS BINARY?
    double E = 0;
    double ESum = 0;
    double ESumSquared = 0;

    // TEMPORARY FOR STEEPEST DESCENT
    double dPsiAlpha = 0;
    double dPsiBeta = 0;
    double dPsiBetaSum = 0;
    double dPsiAlphaSum = 0;
    double dPsiEAlphaSum = 0;
    double dPsiEBetaSum = 0;

    void updateParticle(double **rOld, double **rNew, double &oldWF, double &newWF, int i);

    // USE THIS
    void runMetropolisStep(double **rOld, double **rNew, double &oldWF, double &newWF);
    void sampleSystem(double **r, double newWF, double oldWF);

    WaveFunctions *WF = nullptr;
    MetropolisSampler *R = nullptr;
    void checkRatio(MetropolisSampler *newRatio);
public:
    VMC();
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
