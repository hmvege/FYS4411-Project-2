#ifndef UNIFORMSAMPLING_H
#define UNIFORMSAMPLING_H

#include "metropolissampler.h"

class UniformSampling : public MetropolisSampler
{
private:
    double stepLength;
    std::uniform_real_distribution<double> uniform_distribution;
    double Ratio(double **rPosNew, double **rPosOld, int i, double newWF, double oldWF);
public:
    UniformSampling();
    bool move(double **rPosNew, double **rPosOld, int i, double newWF, double oldWF);
    double nextStep(double **rPosOld, int i, int j);
    void initialize(double newStepLength, double newSeed);
    double initializePosition();
    // Public setters
    void setStepLength(double newStepLength) { stepLength = newStepLength; }
};

#endif // UNIFORMSAMPLING_H
