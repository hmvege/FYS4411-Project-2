#ifndef MetropolisSampler_H
#define MetropolisSampler_H

#include <random>

class MetropolisSampler
{
protected:
    double seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> acceptance_dist; // For choosing to accept a new Metropolis move or not
public:
    MetropolisSampler();
    virtual ~MetropolisSampler();
    virtual bool move(double **rPosNew, double **rPosOld, int i, double newWF, double oldWF);
    virtual double Ratio(double ** rPosNew, double ** rPosOld, int i, double newWF, double oldWF);
    virtual double nextStep(double **rPosOld, int i, int j);
    virtual double initializePosition();
    // Public setters
    void setSeed(double newSeed) { seed = newSeed; }
};

#endif // MetropolisSampler_H
