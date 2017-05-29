#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "metropolissampler.h"
#include "../wavefunctions/wavefunctions.h"

class ImportanceSampler : public MetropolisSampler
{
private:
    double D;
    double deltat;
    double sqrtDeltat;
    double deltatD;
    double **FOld; // Old quantum force
    double **FNew; // New quantum force
    std::normal_distribution<double> gaussian_dist;
//    WaveFunctions *WF = nullptr; // Will this create a double instance of the wavefunction, as we have one stored in the vmc?
    double Ratio(double **rOld, double **rNew, int i, double newWF, double oldWF);
public:
    ImportanceSampler(int new_nParticles, int new_nDimensions, WaveFunctions *newWF);
    ~ImportanceSampler();
    bool move(double **rOld, double **rNew, int i, double newWF, double oldWF);
    void initializePositions(double **rOld, double **rNew);
    void updatePositions(double ** rOld, double ** rNew, int k);
    void initializeSampling(double newStepLength, double newSeed, double newD);
    double GreensRatio(double **y, double **x, int k);
    void setStepLength(double newStepLength) { deltat = newStepLength; }

    // TEMP
    void printQMForces();
};

#endif // IMPORTANCESAMPLER_H
