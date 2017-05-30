#ifndef MetropolisSampler_H
#define MetropolisSampler_H

#include <random>
#include "wavefunctions/wavefunctions.h"

class MetropolisSampler
{
protected:
    int nParticles;
    int nDimensions;
    double seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> acceptance_dist; // For choosing to accept a new Metropolis move or not
    WaveFunctions *WF = nullptr; // Will this create a double instance of the wavefunction, as we have one stored in the vmc?
public:
    MetropolisSampler(int new_nParticles, int new_nDimensions, WaveFunctions *newWF);
    virtual ~MetropolisSampler();
    virtual bool move(double **rOld, double **rNew, int i, double newWF, double oldWF);
    virtual double Ratio(double ** rOld, double ** rNew, int i, double newWF, double oldWF);
    virtual void updatePositions(double ** rOld, double ** rNew, int k);
    virtual void initializePositions(double **rOld, double **rNew);
    // Public setters
    void setSeed(double newSeed) { seed = newSeed; }
    void setWaveFunction(WaveFunctions *newWF) { WF = newWF; }

    // TEMP
    virtual void printQMForces();
};

#endif // MetropolisSampler_H
