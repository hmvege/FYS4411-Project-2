#include "importancesampler.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

ImportanceSampler::ImportanceSampler(int new_nParticles, int new_nDimensions, WaveFunctions *newWF) : MetropolisSampler(new_nParticles, new_nDimensions, newWF)
{
    // Allocating memory to the force matrices.
    FOld = new double * [nParticles];
    FNew = new double * [nParticles];
}

ImportanceSampler::~ImportanceSampler()
{

}

void ImportanceSampler::initializeSampling(double newStepLength, double newSeed, double newD)
{
    /*
     * Initialization of the importance sampler.
     * Arguments:
     * newStepLength    : The step length of the Metropolis update
     * newSeed          : The seed to be used by the Mersenne twister RNG
     * newD             : 0.5 in atomic units
     */
    // Initializing parameters and often used constants
    setStepLength(newStepLength);
    setSeed(newSeed);
    D                   = newD;
    deltatD             = D*deltat;
    sqrtDeltat          = sqrt(deltat);
    // Initializing random generators and distributions
    std::mt19937_64 gen(seed); // Mersenne-Twister19937
    std::uniform_real_distribution <double> uni_dist(0.0,1.0);
    std::normal_distribution <double> gauss_dist(0.0,1.0);
    generator           = gen;
    gaussian_dist       = gauss_dist;
    acceptance_dist     = uni_dist;
}

void ImportanceSampler::initializePositions(double **rOld, double **rNew)
{
    /*
     * Class instance for initializing the electron positions.
     * Arguements:
     * rOld    : Old positions of the electrons
     * rNew    : New positions of the electrons
     */
    for (int i = 0; i < nParticles; i++)
    {
        rOld[i] = new double [nDimensions];
        rNew[i] = new double [nDimensions];
        FOld[i] = new double [nDimensions];
        FNew[i] = new double [nDimensions];
        for (int j = 0; j < nDimensions; j ++)
        {
            rOld[i][j] = gaussian_dist(generator)*sqrtDeltat;
            rNew[i][j] = rOld[i][j];
        }
    }
    WF->initializeWFSampling(rOld);
    for (int i = 0; i < nParticles; i++)
    {
        WF->quantumForce(rOld,FOld,i);
        FNew[i] = FOld[i];
    }
}

void ImportanceSampler::updatePositions(double **rOld, double **rNew, int k)
{
    /*
     * Class instance for updating the position of a single electron.
     * rOld    : Old positions of the electrons
     * rNew    : New positions of the electrons to be determined
     * k       : Particle to be updated
     */
    for (int i = 0; i < nDimensions; i++)
    {
        rNew[k][i] = rOld[k][i] + deltatD*FOld[k][i] + sqrtDeltat*gaussian_dist(generator);
    }
}

bool ImportanceSampler::move(double **rOld, double **rNew, int i, double newWF, double oldWF)
{
    /*
     * Comparing if move should be accepted against an uniform distribution. If accepted, stores the new quantum force as the old.
     * Arguments:
     *  rOld        : old particle positions
     *  rNew        : updated particle positions
     *  i           : particle being moved
     *  newWF       : updated wave function
     *  oldWF       : old wave function
     */
    if (acceptance_dist(generator) <= Ratio(rOld, rNew, i, newWF, oldWF))
    {
        FOld = FNew;
        return true;
    }
    else
    {
        return false;
    }
}

double ImportanceSampler::Ratio(double **rOld, double **rNew, int i, double newWF, double oldWF)
{
    /*
     * Calculates the ratio used for comparing of a move should be accepted or not.
     * Arguments:
     *  rOld        : old particle positions
     *  rNew        : updated particle positions
     *  i           : particle being moved
     *  newWF       : updated wave function
     *  oldWF       : old wave function
     */
    return GreensRatio(rNew, rOld, i)*newWF*newWF/(oldWF*oldWF); // Put the G's together!
}

double ImportanceSampler::GreensRatio(double **y, double **x, int k)
{
    /*
     * Importance sampling ratio from two Greens functions, q = G(x,y)/G(y,x)
     * y : new positions
     * x : old positions
     * k : particle number
     */
    WF->quantumForce(y,FNew,k); // FOld[0] = Fx(x), FOld[1] = Fy[x]
    double GreensFunction = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        GreensFunction += 0.5*(FOld[k][j] + FNew[k][j])*(deltatD*0.5*(FOld[k][j] - FNew[k][j]) - y[k][j] + x[k][j]);
    }
    GreensFunction = exp(GreensFunction);
    return GreensFunction;
}

void ImportanceSampler::printQMForces()
{
    /*
     * Function used for diagnostic and bug-testing purposes.
     */
    cout << "Old QMF" << endl;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            cout << FOld[i][j] << " ";
        }
        cout << endl;
    }
    cout << "New QMF" << endl;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            cout << FNew[i][j] << " ";
        }
        cout << endl;
    }
}
