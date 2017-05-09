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
    return distribution(generator) <= Ratio(rPosNew, rPosOld, i, newWF, oldWF);
}

double UniformSampling::nextStep(double ** rPosOld, int i, int j)
{
    return stepLength * distribution(generator) - 0.5;
}

void UniformSampling::initialize(double newStepLength, double newSeed)
{
    setSeed(newSeed);
    setStepLength(newStepLength);
    std::mt19937_64 gen(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uni_dist(0,1);
    generator = gen;
    uniform_dist = uni_dist;
}

double UniformSampling::initializePosition()
{
    return stepLength * distribution(generator) - 0.5;
}
