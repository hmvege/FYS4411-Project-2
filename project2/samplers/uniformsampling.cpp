#include "uniformsampling.h"
#include "wavefunctions/wavefunctions.h"

#include "iostream"

UniformSampling::UniformSampling(int new_nParticles, int new_nDimensions) : MetropolisSampler(new_nParticles, new_nDimensions)
{

}

bool UniformSampling::move(double ** rOld, double ** rNew, int i, double newWF, double oldWF)
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
        return true;
    }
    else
    {
        return false;
    }
}

double UniformSampling::Ratio(double ** rOld, double ** rNew, int i, double newWF, double oldWF)
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
    return (newWF*newWF)/(oldWF*oldWF);
}

void UniformSampling::updatePositions(double **rOld, double **rNew, int k)
{
    /*
     * Class instance for updating the position of a single electron.
     * rOld    : Old positions of the electrons
     * rNew    : New positions of the electrons to be determined
     * k       : Particle to be updated
     */
    for (int i = 0; i < nDimensions; i++)
    {
        rNew[k][i] = rOld[k][i] + stepLength*uniform_distribution(generator);
    }
}

void UniformSampling::initializeSampling(double newStepLength, double newSeed)
{
    /*
     * Initialization of the uniform sampler.
     * Arguments:
     * newStepLength    : The step length of the Metropolis update
     * newSeed          : The seed to be used by the Mersenne twister RNG
     */
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
    /*
     * Class instance for initializing the electron positions.
     * Arguements:
     * rOld    : Old positions of the electrons
     * rNew    : New positions of the electrons
     */
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j ++)
        {
            rOld[i][j] = acceptance_dist(generator)*2.0 - 1.0; // STORE OLD QM FORCE
            rNew[i][j] = rOld[i][j];
        }
    }
    WF->initializeWFSampling(rOld);
}
