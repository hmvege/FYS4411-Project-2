#include "nelectron.h"
#include <iostream>
#include "../wavefunctions.h"
#include "hermite.h"
#include "state.h"
#include "functions.h"

#include <armadillo>

using std::cout;
using std::endl;

// PUT INTO FUNCTION-CONTAINER?
//int factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }

NElectron::NElectron(int new_nParticles, int new_nDimensions, int new_nVarParams, double new_omega, double new_alpha, double new_beta) : WaveFunctions(new_nParticles, new_nDimensions, new_nVarParams)
{
    /*
     * Class for the wavefunction of the N-electron case. Takes only closed shell systems.
     * Arguments:
     * -------------------------------------------------------------------------
     * TODO:
     * [x] implement a usable settup/mapping for the different electron states.
     * [x] implement wavefunctions
     * [x] initialize the set up of the slater determinants
     * [ ] implement matrix inverse
     * [ ] implement inverse update
     * [ ] implement wf functions
     * -------------------------------------------------------------------------
     */
    setOmega(new_omega);
    setAlpha(new_alpha);
    setBeta(new_beta);
    // Hardcoded the maxshell -> BAD
    int maxShell;
    if (nParticles==2) maxShell = 1;
    else if (nParticles==6) maxShell = 2;
    else if (nParticles==12) maxShell = 3;
    else
    {
        printf("Error: %d does not form a closed shell system of electron.\n", nParticles);
        exit(1);
    }
    // Initializing states mappin
    int j_states;
    states = new State[nParticles];
    j_states = 0; // for quantum states
    for (int shellInt = 0; shellInt < maxShell; shellInt++) // Ensuring we are adding the right shell
    {
        for (int n_x = 0; n_x < maxShell; n_x++) // Running over variations of nx
        {
            for (int n_y = 0; n_y < maxShell; n_y++) // Running over variations of ny
            {
                for (double spin = -1; spin < 2; spin += 2) // Running over the two possible spin configurations
                {
                    if (n_x+n_y == shellInt) // If-test to prevent counting outside of the shell
                    {
                        if (spin==-1)
                        {
                            states[j_states].set(n_x,n_y,spin);
                            states[j_states].setHermite(&hermite);
                        }
                        else
                        {
                            states[j_states].set(n_x,n_y,spin);
                            states[j_states].setHermite(&hermite);
                        }
                        j_states++;
                    }
                }
            }
        }
    }

//    // Testing
//    for (int i = 0; i < nParticles; i++)
//    {
//        cout << "i=" << i << ", State: ";
//        states[i].print();
//    }
}

NElectron::~NElectron()
{
    // FIX DESTRUCTOR
}

void NElectron::initialize(double **r)
{
    /*
     * Initializes the slater determinant for our wave function.
     * Arguments:
     *  r   : particle positions
     */
    initializeSlater(r);
}

double NElectron::calculate(double **r)
{
    /*
     * Returns the wave function for N electrons.
     * Arguments:
     *  r   : particle positions
     */
    return psiJastrow(r)*psiSlater(r);
}


double NElectron::localEnergy(double **r)
{
    /*
     * Returns the local energy for N electrons.
     * Arguments:
     *  r   : particle positions
     */
    double energy = 0;
    for (int i = 0; i < nParticles; i++)
    {
        energy += -0.5*laplacian(r,i) + 0.5*omega*omega*(r[i][0] + r[i][1])*(r[i][0] + r[i][1]);
    }
    if (coulombInteraction)
    {
        return energy + coulomb(r);
    }
    else
    {
        return energy;
    }
}

void NElectron::quantumForce(double **r, double **F, int k)
{
    /*
     * Returns the quantum force for N electrons.
     * Arguments:
     *  r   : particle positions
     */
    //Update slater determinants here?

}

void NElectron::steepestDescent(double &E, int NCycles)
{

}

void NElectron::sampleSD(double **r, double &E)
{

}

bool NElectron::SDConvergenceCriteria()
{
    return false;
}

void NElectron::initializeSlater(double **r)
{
    /*
     * Sets the a spin up and spin down Slater matrix
     */
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            DSpinDown[i][j] = states[j].wf(r[i], alpha, omega);
            DSpinUp[i][j] = states[j*2].wf(r[i*2], alpha, omega);
        }
    }

}

void NElectron::updateSlater(double **r, int k)
{
    /*
     * Updates the inverse of the slater determinant as well as the slater determinant itself.
     * Arguments:
     *  r   : particle positions
     *  k   : particle being moved
     */
}

//double NElectron::determinant(double **A, int dim)
//{
//    /*
//     * Function for calculating the determinant.
//     * Taken from lecture.
//     * Arguments:
//     *  A   : square matrix
//     *  dim : dimensionality of the matrix
//     */
//    if (dim == 2)
//    {
//        return A[0][0]*A[1][1] - A[0][1]*A[1][0];
//    }
//    double sum = 0;
//    for (int i = 0; i < dim; i++)
//    {
//        double ** sub = new double*[dim-1];
//        for (int j = 0; j < i; j++)
//        {
//            sub[j] = &A[j][1];
//        }
//        for (int j = i + 1; j < dim; j++)
//        {
//            sub[j-1] = &A[j][1];
//        }
//        if (i % 2 == 0)
//        {
//            sum += A[i][0] * determinant(sub, dim-1);
//        }
//        else
//        {
//            sum -= A[i][0] * determinant(sub, dim-1);
//        }
//    }
//    return sum;
//}

double NElectron::psiJastrow(double **r)
{
    /*
     * Returns the correlation term(the Jastrow factor) for the particles.
     * Arguments:
     *  r   : position
     */
    double psiSum = 0;
    for (int i = 0; i < nParticles; i++) // running over particles
    {
        for (int j = 0; j < i; j++)
        {
            if (states[i].getSpin() + states[j].getSpin() == 0) // a = 1.0
            {
                psiSum -= 1.0 / (1.0/r_ij(r[i],r[j]) + beta);
            }
        }
    }
    return exp(psiSum);
}

double NElectron::psiSlater(double **r)
{
    /*
     * Returns the Slater determninant for the electrons.
     * Arguments:
     *  r   : position
     */
    return determinant(DSpinUp, nParticles)*determinant(DSpinDown, nParticles)/factorial(nParticles/2); // SHOULD IT BE HALF HERE OR sqrt(N!)?
}

void NElectron::gradientJastrow(double * grad, double **r, int k)
{
    /*
     * Finds the gradient of the Jastrow wave function.
     * Arguments:
     *  grad    : reference to gradient array of length 2
     *  r       : positions of the particles
     *  k       : particle we are getting the gradient for
     */
    grad[0] = 0;
    grad[1] = 0;
    double a = 0;
    double commonFactor = 0;
    double r_dist = 0;
    for (int i = 0; i < k; i++) // STOP AT k-1 OR k?
    {
        a = checkSpin(i,k);
        r_dist = r_ij(r[i],r[k]);
        commonFactor = a/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
        grad[0] += (r[i][0]-r[k][0])*commonFactor; // Is the difference r_i and r_k correct?
        grad[1] += (r[i][1]-r[k][1])*commonFactor;
    }
    for (int i = k+1; i < nParticles; i++)
    {
        a = checkSpin(k,i); //
        r_dist = r_ij(r[k],r[i]);
        commonFactor = a/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
        grad[0] += (r[k][0]-r[i][0])*commonFactor; // Is the difference r_i and r_k correct?
        grad[1] += (r[k][1]-r[i][1])*commonFactor;
    }
}

void NElectron::gradientSD(double * grad, double **r, int k)
{
    /*
     * Finds the gradient of the Slater determinant.
     * Arguments:
     *  grad    : reference to gradient array of length 2
     *  r       : positions of the particles
     *  k       : particle we are getting the gradient for
     */
    double *wfGrad = new double[2];
    wfGrad[0] = 0;
    wfGrad[1] = 0;
    for (int i = 0; i < nParticles; i++)
    {
        states[i].wfGradient(wfGrad,r[k],alpha,omega);
        grad[0] += wfGrad[0]*SpinUpInverse[i][k];
        grad[1] += wfGrad[1]*SpinDownInverse[i][k]; // ??????????????????????????????????????????????????????????????????????????????
    }
    delete [] wfGrad;
}

double NElectron::laplacianJastrow(double **r, int k)
{
    /*
     * Returns the laplacian of the Jastrow factor.
     * Arguments:
     *  r   : particle positions
     *  k   : index of particle being moved
     */
    double lap = 0;
    double a = 0;
    double r_dist = 0;
    // Method 1
    double sum1x = 0;
    double sum1y = 0;
    double sum2x = 0;
    double sum2y = 0;
    double commonFactor = 0;
    double r_ijBeta = 0;
    for (int i = 0; i < k; i++)  // STOP AT k-1 OR k?
    {
        a = checkSpin(i,k);
        r_dist = r_ij(r[i],r[k]);
        r_ijBeta = (1+beta*r_dist);
        commonFactor = a/(r_dist*r_ijBeta*r_ijBeta);
        sum1x += (r[i][0]-r[k][0])*commonFactor; // Is the difference r_i and r_k correct?
        sum1y += (r[i][1]-r[k][1])*commonFactor;
    }
    for (int i = k+1; i < nParticles; i++)
    {
        a = checkSpin(k,i);
        r_dist = r_ij(r[k],r[i]);
        r_ijBeta = (1+beta*r_dist);
        commonFactor = a/(r_dist*r_ijBeta*r_ijBeta);
        sum2x += (r[k][0]-r[i][0])*commonFactor; // Is the difference r_i and r_k correct?
        sum2y += (r[k][1]-r[i][1])*commonFactor;
    }
    lap += (sum1x+sum2x)*(sum1x+sum2x) + (sum1y+sum2y)*(sum1y+sum2y);
    for (int i = 0; i < k; i++)
    {
        a = checkSpin(i,k);
        r_dist = r_ij(r[i],r[k]);
        r_ijBeta = (1+beta*r_dist);
        commonFactor = a/(r_ijBeta*r_ijBeta);
        lap += 1.0/r_dist*commonFactor - 2*beta*commonFactor/r_ijBeta;
    }
    for (int i = k+1; i < nParticles; i++)  // STOP AT k-1 OR k?
    {
        a = checkSpin(k,i);
        r_dist = r_ij(r[k],r[i]);
        r_ijBeta = (1+beta*r_dist);
        commonFactor = a/(r_ijBeta*r_ijBeta);
        lap += 1.0/r_dist*commonFactor - 2*beta*commonFactor/r_ijBeta;
    }
    // Method 2
//    double commonFactor_i = 0;
//    double commonFactor_j = 0;
//    for (int i = 0; i < nParticles; i++)
//    {
//        if (i==k) continue;
//        a = checkSpin(k,i);
//        r_dist = r_ij(r[k],r[i]);
//        commonFactor_i = a/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
//        for (int j = 0; j < nParticles; j++)
//        {
//            if (j==k) continue;
//            r_dist = r_ij(r[k],r[j]);
//            commonFactor_j = checkSpin(k,j)/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
//            lap += ((r[k][0] - r[i][0])*(r[k][0] - r[j][0]) + (r[k][1] - r[i][1])*(r[k][1] - r[j][1]))*commonFactor_i*commonFactor_j;
//        }
//    }
//    for (int i = 0; i < nParticles; i++)
//    {
//        if (i==k) continue;
//        r_dist = r_ij(r[k],r[i]);
//        commonFactor_i = 1 + beta*r_dist;
//        lap += 2*checkSpin(k,i)/(r_dist*commonFactor_i*commonFactor_i*commonFactor_i);
//    }
    return lap;
}

double NElectron::laplacianSD(double **r, int k)
{
    double lap = 0;
    for (int i = 0; i < nParticles; i++)
    {
        lap += states[i].wfLaplacian(r[k], alpha, omega) * DInv[i][k];
    }
    return lap;
}

double NElectron::laplacian(double **r, int k)
{
    /*
     * Returns the laplacian of a particle.
     * Arguments:
     *  r   : particle positions
     *  k   : particle to find the laplacian for
     */
    double *gradSlater = new double[2];
    double *gradJastrow = new double[2];
    gradSlater[0] = 0;
    gradSlater[1] = 0;
    gradJastrow[0] = 0;
    gradJastrow[1] = 0;

    // ??????????????????????????????????????????????????????????????????????????????

    delete [] gradSlater;
    delete [] gradJastrow;

    return 0.0;
}

double NElectron::checkSpin(int i, int j)
{
    /*
     * Checks the spin of between two particles.
     * Arguments:
     *  i    : particle i
     *  k    : particle k
     * Returns:
     *  1.0 if spins are anti-parallel
     *  1/3 if spins are parallel
     */
    if (states[i].getSpin() + states[j].getSpin() == 0) // Returns 1 if spins are anti-parallel
    {
        return 1.0;
    }
    else // Returns 1/3 if spins are parallel
    {
        return 0.333333333333333;
    }
}
