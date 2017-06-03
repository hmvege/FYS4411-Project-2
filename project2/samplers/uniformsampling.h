#ifndef UNIFORMSAMPLING_H
#define UNIFORMSAMPLING_H

#include "metropolissampler.h"

class UniformSampling : public MetropolisSampler
{
private:
    double stepLength;
    std::uniform_real_distribution<double> uniform_distribution;
    double Ratio(double **rOld, double **rNew, int i, double newWF, double oldWF);
public:
    UniformSampling(int new_nParticles, int new_nDimensions, int new_numprocs, int new_processRank);
    bool move(double **rOld, double **rNew, int i, double newWF, double oldWF);    
    void initializeSampling(double newStepLength, double newSeed);
    void initializePositions(double **rOld, double **rNew);
    void updatePositions(double **rOld, double **rNew, int k);
    // Public setters
    void setStepLength(double newStepLength) { stepLength = newStepLength; }
};

#endif // UNIFORMSAMPLING_H
