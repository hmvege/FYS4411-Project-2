#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "metropolissampler.h"
#include "../wavefunctions/wavefunctions.h"

class ImportanceSampler : public MetropolisSampler
{
private:
    int nParticles;
    int nDimensions;
    double D;
    double deltat;
    double sqrtDeltat;
    double deltatD;
    double exp_denom_factor;
    double denom_factor;
    // ADD OLD QM FORCE
    // ADD NEW QM FORCE
    std::normal_distribution<double> gaussian_dist;
    WaveFunctions *WF = nullptr; // Will this create a double instance of the wavefunction, as we have one stored in the vmc?
    double Ratio(double **rPosNew, double **rPosOld, int i, double newWF, double oldWF);
public:
    ImportanceSampler();
    ~ImportanceSampler();
    bool move(double **rPosNew, double **rPosOld, int i, double newWF, double oldWF);
    double initializePosition();
    double nextStep(double **rPosOld, int i, int j);
    void initializeSampling(double newStepLength, double newSeed, double newD, int newNPart, int newNDim);
    double q(double **y, double **x, int k);
    void setWaveFunction(WaveFunctions *newWF) { WF = newWF; }
    void setStepLength(double newStepLength) { deltat = newStepLength; }
};

#endif // IMPORTANCESAMPLER_H
