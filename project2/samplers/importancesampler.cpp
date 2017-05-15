#include "importancesampler.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

ImportanceSampler::ImportanceSampler(int new_nParticles, int new_nDimensions) : MetropolisSampler(new_nParticles, new_nDimensions)
{
    // Allocating memory to the force matrices.
    FOld = new double * [nParticles];
    FNew = new double * [nParticles];
}

ImportanceSampler::~ImportanceSampler()
{

}

double ImportanceSampler::Ratio(double **rOld, double **rNew, int i, double newWF, double oldWF)
{
    return q(rNew, rOld, i)*newWF*newWF/(oldWF*oldWF); // Put the G's together!
}

bool ImportanceSampler::move(double **rOld, double **rNew, int i, double newWF, double oldWF)
{
    /*
     * Comparing if move should be accepted against an uniform distribution. If accepted, stores the new quantum force as the old.
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
//    return acceptance_dist(generator) <= Ratio(rPosNew, rPosOld, i, newWF, oldWF); // CORRECT COMPARISON
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

void ImportanceSampler::initializePositions(double **rOld, double **rNew)
{
    /*
     * Class instance for initializing the electron positions.
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
            rOld[i][j] = acceptance_dist(generator)*2.0 - 1.0; // STORE OLD QM FORCE
            rNew[i][j] = rOld[i][j];
        }
    }
    for (int i = 0; i < nParticles; i++)
    {
        WF->quantumForce(rOld,FOld,i);
        FNew[i] = FOld[i];
    }
//    return acceptance_dist(generator)*2.0 - 1.0; // Ensures we are working on a uniform interaval from -1 to 1
}

void ImportanceSampler::initializeSampling(double newStepLength, double newSeed, double newD)
{
    // Initializing parameters and often used constants
    setStepLength(newStepLength);
    setSeed(newSeed);
    D                   = newD;
    deltatD             = D*deltat;
    sqrtDeltat          = sqrt(deltat);
    // Initializing random generators and distributions
    std::mt19937_64 gen(seed); // Mersenne-Twister19937
    std::uniform_real_distribution <double> uni_dist(0.0,1.0);
    std::normal_distribution <double> gauss_dist(0.0,1.0); // Mean should be around 0 and the std should be 1?
    generator           = gen;
    gaussian_dist       = gauss_dist;
    acceptance_dist     = uni_dist;
}

double ImportanceSampler::q(double **y, double **x, int k)
{
    /*
     * Importance sampling ratio from two Greens functions, q = G(x,y)/G(y,x)
     * y : new positions
     * x : old positions
     * k : particle number
     */
    WF->quantumForce(y,FNew,k); // FOld[0] = Fx(x), FOld[1] = Fy[x]
    return exp(0.5*((x[k][0]-y[k][0])*(FOld[k][0]+FNew[k][0])  +
                    (x[k][1]-y[k][1])*(FOld[k][1]+FNew[k][1])) +
            deltatD*0.25*(FOld[k][0]*FOld[k][0] + FOld[k][1]*FOld[k][1] - FNew[k][1]*FNew[k][1] - FNew[k][0]*FNew[k][0]));
}

// REMOVE ONCE DONE
void ImportanceSampler::printQMForces()
{
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
