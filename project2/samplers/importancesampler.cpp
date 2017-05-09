#include "importancesampler.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

ImportanceSampler::ImportanceSampler()
{

}

ImportanceSampler::~ImportanceSampler()
{

}

double ImportanceSampler::Ratio(double ** rPosNew, double ** rPosOld, int i, double newWF, double oldWF)
{
    return q(rPosNew, rPosOld, i)*newWF*newWF/(oldWF*oldWF); // Put the G's together!
}

bool ImportanceSampler::move(double **rPosNew, double **rPosOld, int i, double newWF, double oldWF)
{
    return uniform_dist(generator) <= Ratio(rPosNew, rPosOld, i, newWF, oldWF);
}

double ImportanceSampler::nextStep(double **rPosOld, int i, int j)
{
    double *F = WF->quantumForce(rPosOld,i); // IS THIS GENERATE MEMORY LEAKAGE?
    return deltatD*F[j] + sqrtDeltat*gaussian_dist(generator);
}

double ImportanceSampler::initializePosition()
{
    return deltat * uniform_dist(generator) - 0.5;
}

void ImportanceSampler::initializeSampling(double newStepLength, double newSeed, double newD, int newNPart, int newNDim)
{
    // Initializing parameters and often used constants
    setStepLength(newStepLength);
    setSeed(newSeed);
    nParticles          = newNPart;
    nDimensions         = newNDim;
    D                   = newD;
    deltatD             = D*deltat;
    sqrtDeltat          = sqrt(deltat);
    exp_denom_factor    = 1.0/(4*deltatD);
    denom_factor        = 1.0/pow(4*M_PI*deltatD, 3*double(nParticles)/2.);
    // Initializing random generators and distributions
    std::mt19937_64 gen(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution <double> uni_dist(0,1.0);
    std::normal_distribution <double> gauss_dist(0,0.5); // Mean should be around 0 and the std should be what? ??????
    cout << "importancesampler.cpp: Remember to set the normal distribution correctly" << endl;
    generator           = gen;
    gaussian_dist       = gauss_dist;
    uniform_dist        = uni_dist;
}

double ImportanceSampler::q(double **y, double **x, int k)
{
    /*
     * Importance sampling ratio from two Greens functions, q = G(x,y)/G(y,x)
     * y : new positions
     * x : old positions
     * k : particle number
     */
    double *FOld = WF->quantumForce(x,k);
    double *FNew = WF->quantumForce(y,k);
    return exp(0.5*(y[k][0]-x[k][0])*(FOld[0]+FNew[0]) + 0.5*(y[k][1]-x[k][1])*(FOld[1]+FNew[1]) + deltat*0.25*(FNew[0]*FNew[0] - FOld[0]*FOld[0] + FNew[1]*FNew[1] - FOld[1]*FOld[1]));
}
