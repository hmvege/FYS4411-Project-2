#include "uniformsampling.h"

UniformSampling::UniformSampling(int new_nParticles, int new_nDimensions) : MetropolisSampler(new_nParticles, new_nDimensions)
{
//    cout << "Remember to comment on uniform sampling" << edl;
}

double UniformSampling::Ratio(double ** rOld, double ** rNew, int i, double newWF, double oldWF)
{
    return (newWF*newWF)/(oldWF*oldWF);
}

bool UniformSampling::move(double ** rOld, double ** rNew, int i, double newWF, double oldWF)
{
    return acceptance_dist(generator) <= Ratio(rOld, rNew, i, newWF, oldWF);
}

double UniformSampling::nextStep(double ** rOld, int i, int j)
{
    return  stepLength * uniform_distribution(generator);
}

void UniformSampling::updatePositions(double **rOld, double **rNew, int k)
{
    for (int i = 0; i < nDimensions; i++)
    {
        rNew[k][i] = rOld[k][i] + stepLength*uniform_distribution(generator);
    }
}

void UniformSampling::initialize(double newStepLength, double newSeed)
{
    setSeed(newSeed);
    setStepLength(newStepLength);
    std::mt19937_64 gen(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uni_dist(-1,1);
    std::uniform_real_distribution<double> accept_dist(0,1);
    generator = gen;
    uniform_distribution = uni_dist;
    acceptance_dist = accept_dist;
}

void UniformSampling::initializePositions(double **rOld, double **rNew)
{
    for (int i = 0; i < nParticles; i++)
    {
        rOld[i] = new double[nDimensions];
        rNew[i] = new double[nDimensions];
        for (int j = 0; j < nDimensions; j ++)
        {
            rOld[i][j] = acceptance_dist(generator)*2.0 - 1.0; // STORE OLD QM FORCE
            rNew[i][j] = rOld[i][j];
        }
    }
//    return uniform_distribution(generator);
}
