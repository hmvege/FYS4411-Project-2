#include "uniformsampling.h"

UniformSampling::UniformSampling()
{

}

double UniformSampling::Ratio(double ** rPosNew, double ** rPosOld, int i, double newWF, double oldWF)
{
    return (newWF*newWF)/(oldWF*oldWF);
}

bool UniformSampling::move(double ** rPosNew, double ** rPosOld, int i, double newWF, double oldWF)
{
    return acceptance_dist(generator) <= Ratio(rPosNew, rPosOld, i, newWF, oldWF);
}

double UniformSampling::nextStep(double ** rPosOld, int i, int j)
{
    return stepLength * uniform_distribution(generator);
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

double UniformSampling::initializePosition()
{
    return uniform_distribution(generator);
}
