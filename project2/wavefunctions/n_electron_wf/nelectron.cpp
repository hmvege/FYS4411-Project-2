#include "nelectron.h"
#include <iostream>
#include "../wavefunctions.h"
#include "hermite.h"
#include "state.h"
#include "functions.h"

#include <armadillo>

// TEMP
#include <iomanip>

using std::cout;
using std::endl;

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
     * [x] implement matrix inverse
     * [ ] implement inverse update
     * [x] implement wf functions
     * [ ] implement quantum force
     * [ ] implement steepest descent
     * [x] get correct results for 2 electrons->identify and remove bugs!
     * -------------------------------------------------------------------------
     */
    setOmega(new_omega);
    setAlpha(new_alpha);
    setBeta(new_beta);
    // Hardcoded the maxshell input -> BAD? Try to find a mapping of some sort
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
    states = new State*[nParticles];
    for (int i = 0; i < nParticles; i++)
    {
        states[i] = new State;
    }
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
                            states[j_states]->set(n_x,n_y,spin);
                            states[j_states]->setHermite(&hermite);
                        }
                        else
                        {
                            states[j_states]->set(n_x,n_y,spin);
                            states[j_states]->setHermite(&hermite);
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
//        states[i]->print();
//    }
}

NElectron::~NElectron()
{
    for (int i = 0; i < nParticles/2; i++)
    {
        delete [] DSpinDown[i];
        delete [] DSpinUp[i];
        delete [] DSpinDownInverse[i];
        delete [] DSpinUpInverse[i];
        delete [] DSpinDownOld[i];
        delete [] DSpinUpOld[i];
        delete [] DSpinDownInverseOld[i];
        delete [] DSpinUpInverseOld[i];
    }
    delete [] DSpinDown;
    delete [] DSpinUp;
    delete [] DSpinDownInverse;
    delete [] DSpinUpInverse;
    delete [] DSpinDownOld;
    delete [] DSpinUpOld;
    delete [] DSpinDownInverseOld;
    delete [] DSpinUpInverseOld;
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
    updateSlater(r); // UPDATE SLATER HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // TEMPORARY
    double slaterWF = psiSlater(r);
    if (runJastrow)
    {
        double jastrowWF = psiJastrow(r);
        return jastrowWF*slaterWF;
    }
    else
    {
        return slaterWF;
    }
//    return psiJastrow(r)*psiSlater(r);
}


double NElectron::localEnergy(double **r)
{
    /*
     * Returns the local energy for N electrons.
     * Arguments:
     *  r   : particle positions
     */
//    updateSlater(r);
    double energy = 0;
    for (int i = 0; i < nParticles; i++)
    {
        energy += -0.5*laplacian(r,i) + 0.5*omega*omega*(r[i][0]*r[i][0] + r[i][1]*r[i][1]);
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
     *  F   : quantum forces for all particles
     *  k   : particle we are finding the quantum force for
     */
    double *gradSlater = new double[2];
    double *gradJastrow = new double[2];
    gradSlater[0] = 0;
    gradSlater[1] = 0;
    gradJastrow[0] = 0;
    gradJastrow[1] = 0;
//    printf("k=%2d grad=( %10f , %10f ) r[k] = ( %10f , %10f ) \n", k, gradSlater[0],gradSlater[1],r[k][0],r[k][1]);
    // Running with/without Jastrow
    if (runJastrow)
    {
        gradientJastrow(gradJastrow,r,k);
        gradientSlater(gradSlater,r,k);
        F[k][0] = 2*(gradSlater[0] + gradJastrow[0]);
        F[k][1] = 2*(gradSlater[1] + gradJastrow[1]);
    }
    else
    {
        gradientSlater(gradSlater,r,k);
        F[k][0] = 2*gradSlater[0];
        F[k][1] = 2*gradSlater[1];
    }

    delete [] gradSlater;
    delete [] gradJastrow;
}

void NElectron::steepestDescent(double &ESum, int NCycles)
{
    /*
     * Using steepest descent to update the variational parameters of the wavefunction.
     * Arguments:
     *  ESum        : Local energy sum from running N MC cycles
     *  NCyclces    : Number of MC cycles used
     */
    SDStatistics(NCycles);
    double alphaDerivative = 2*(dPsiEAlphaSum - dPsiAlphaSum*ESum);
    double betaDerivative = 2*(dPsiEBetaSum - dPsiBetaSum*ESum);
    alpha -= SDStepLength*alphaDerivative; // Updating alpha and beta
    beta -= SDStepLength*betaDerivative;
}

void NElectron::sampleSD(double **r, double &E)
{
    /*
     * Sampling used by the steepest descent algorithm.
     * Arguments:
     *  r   : particle positions
     *  E   : local energy of current positions
     */
    //// PART TO CHANGE =============================================================================================================
    dPsiAlpha       = 0; // Derivative of WF w.r.t. alpha
    dPsiBeta        = 0; // Derivative of WF w.r.t. beta
    //// ============================================================================================================================
    dPsiAlphaSum    += dPsiAlpha;
    dPsiEAlphaSum   += dPsiAlpha*E;
    dPsiBetaSum     += dPsiBeta;
    dPsiEBetaSum    += dPsiBeta*E;
}

void NElectron::SDStatistics(int NCycles)
{
    /*
     * Function for retrieving Steepest Descent statistics. Arguments:
     * r        : positions
     * NCycles  : Monte Carlo cycles
     */
    dPsiBetaSum     /= double(nParticles*NCycles);
    dPsiAlphaSum    /= double(nParticles*NCycles);
    dPsiEAlphaSum   /= double(nParticles*NCycles);
    dPsiEBetaSum    /= double(nParticles*NCycles);
}

void NElectron::printVariationalParameters()
{
    /*
     * Temporary function for printing the variational parameters used.
     */
    cout << "Alpha = " << std::setw(10) << alpha << " Beta = " << std::setw(10) << beta << endl;
}

bool NElectron::SDConvergenceCriteria()
{
    printf("Steepest Descent convergence criteria not implemented\n");
    return false;
}

void NElectron::initializeSlater(double **r)
{
    /*
     * Sets the a spin up and spin down Slater matrix. Please ignore the terrible hard-coded matrix allocation.
     */
//    DSpinDown = arma::zeros<arma::mat>(nParticles/2,nParticles/2);
//    DSpinUp = arma::zeros<arma::mat>(nParticles/2,nParticles/2);
//    DSpinDownInverse = arma::zeros<arma::mat>(nParticles/2,nParticles/2);
//    DSpinUpInverse = arma::zeros<arma::mat>(nParticles/2,nParticles/2);
    DSpinDown           = new double*[nParticles/2];
    DSpinUp             = new double*[nParticles/2];
    DSpinDownInverse    = new double*[nParticles/2];
    DSpinUpInverse      = new double*[nParticles/2];
    DSpinDownOld        = new double*[nParticles/2];
    DSpinUpOld          = new double*[nParticles/2];
    DSpinDownInverseOld = new double*[nParticles/2];
    DSpinUpInverseOld   = new double*[nParticles/2];
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        DSpinDown[i]            = new double[nParticles/2];
        DSpinUp[i]              = new double[nParticles/2];
        DSpinDownInverse[i]     = new double[nParticles/2];
        DSpinUpInverse[i]       = new double[nParticles/2];
        DSpinDownOld[i]         = new double[nParticles/2];
        DSpinUpOld[i]           = new double[nParticles/2];
        DSpinDownInverseOld[i]  = new double[nParticles/2];
        DSpinUpInverseOld[i]    = new double[nParticles/2];
        for (int j = 0; j < nParticles/2; j++) // States
        {
//            printf("Spin down: (%5d,%5d) Spin up: (%5d%5d)\n", 2*j, 2*i, 2*j+1, 2*i+1);
            DSpinDown[i][j]             = states[2*j]->wf(r[2*i], alpha, omega);
            DSpinUp[i][j]               = states[2*j+1]->wf(r[2*i+1], alpha, omega);
            DSpinDownInverse[i][j]      = DSpinDown[i][j];
            DSpinUpInverse[i][j]        = DSpinUp[i][j];
//            DSpinDownOld[i][j]          = DSpinDown[i][j];
//            DSpinUpOld[i][j]            = DSpinUp[i][j];
//            DSpinDownInverseOld[i][j]   = DSpinDownInverse[i][j];
//            DSpinUpInverseOld[i][j]     = DSpinUpInverse[i][j];
//            DSpinDown(i,j) = states[2*j]->wf(r[2*i], alpha, omega);
//            DSpinUp(i,j) = states[2*j+1]->wf(r[2*i+1], alpha, omega);
//            DSpinDownInverse(i,j) = DSpinDown(i,j);
//            DSpinUpInverse(i,j) = DSpinUp(i,j);
        }
    }
//    DSpinDownInverse.i();
//    DSpinUpInverse.i();
    inverse(DSpinDownInverse, nParticles/2);
    inverse(DSpinUpInverse, nParticles/2);
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            DSpinDownOld[i][j]          = DSpinDown[i][j];
            DSpinUpOld[i][j]            = DSpinUp[i][j];
            DSpinDownInverseOld[i][j]   = DSpinDownInverse[i][j];
            DSpinUpInverseOld[i][j]     = DSpinUpInverse[i][j];
        }
    }
}

void NElectron::updateSlater(double **r)
{
    /*
     * Updates the inverse of the slater determinant as well as the slater determinant itself.
     * Arguments:
     *  r   : particle positions
     *  k   : particle being moved
     */
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            // Storing old slater matrices
            DSpinDownOld[i][j]          = DSpinDown[i][j];
            DSpinUpOld[i][j]            = DSpinUp[i][j];
            DSpinDownInverseOld[i][j]   = DSpinDownInverse[i][j];
            DSpinUpInverseOld[i][j]     = DSpinUpInverse[i][j];
            // Getting new slater matrices
            DSpinDown[i][j]             = states[2*j]->wf(r[2*i], alpha, omega);
            DSpinUp[i][j]               = states[2*j+1]->wf(r[2*i+1], alpha, omega);
            DSpinDownInverse[i][j]      = DSpinDown[i][j];
            DSpinUpInverse[i][j]        = DSpinUp[i][j];
//            DSpinDown(i,j) = states[2*j]->wf(r[2*i], alpha, omega);
//            DSpinUp(i,j) = states[2*j+1]->wf(r[2*i+1], alpha, omega);
//            DSpinDownInverse(i,j) = DSpinDown(i,j);
//            DSpinUpInverse(i,j) = DSpinUp(i,j);
        }
    }
//    DSpinDownInverse.i();
//    DSpinUpInverse.i();
    // Finding the inverse matrices - COULD BE DONE MORE EFFICIENTLY!!
    inverse(DSpinDownInverse, nParticles/2);
    inverse(DSpinUpInverse, nParticles/2);
//    printf("D = %10f D_inverse = %10f D*D_inverse = %10f \n", DSpinDown[0][0], DSpinDownInverse[0][0], DSpinDown[0][0]*DSpinDownInverse[0][0]);
//    if (k%2==0) // Spin down
//    {
//        for (int i = 0; i < nParticles/2; i++) // States
//        {
//            DSpinDown[k][i] = states[2*k].wf(r[2*i], alpha, omega);
//        }
//    }
//    else // Spin up
//    {
//        for (int i = 0; i < nParticles/2; i++) // States
//        {
//            DSpinUp[k][i] = states[2*k+1].wf(r[2*i+1], alpha, omega);
//        }
//    }
    // Updating the inverse of the Slater matrix
//    double R = 0; // Slater determinant ratio
//    for (int i = 0; i < nParticles/2; i++)
//    for (int i = 0; i < nParticles/2; i++)
//    {
//        for (int j = 0; j < nParticles/2; j++)
//        {
//            if (j!=i)
//            {
//                DSpinDownInverse[i][j] -= d[i][k]/R
//            }
//            else
//            {
//                DSpinDownInverse[i][j] = d[i][k]/R *
//            }
//        }
//    }
}

double NElectron::updateInverseSlater(double **r, int k)
{
    /*
     * An efficient way of updating the inverse of the Slater matrices.
     */
//    double sum = 0;
//    for (int l = 0; l < nParticles; l++)
//    {

//    }
    return 0.0;
}

double NElectron::psiSlater(double **r)
{
    /*
     * Returns the Slater determninant for the electrons.
     * Arguments:
     *  r   : position
     */
//    return determinant(DSpinUp, nParticles/2)*determinant(DSpinDown, nParticles/2)/double(factorial(nParticles/2)); // SHOULD IT BE HALF HERE OR sqrt(N!)?
    return determinant(DSpinUp, nParticles/2)*determinant(DSpinDown, nParticles/2); // Since we divide wavefunctions on each other, we do not need factorial
//    return det(DSpinUp, nParticles/2)*det(DSpinDown, nParticles/2); // Since we divide wavefunctions on each other, we do not need factorial
}

void NElectron::gradientSlater(double * grad, double **r, int k)
{
    /*
     * Finds the gradient of the Slater determinant.
     * Arguments:
     *  grad    : reference to gradient array of length 2
     *  r       : positions of the particles
     *  k       : particle we are getting the gradient for
     */
//    double *wfGradUp = new double[2];
//    double *wfGradDown = new double[2];
//    wfGradUp[0] = 0;
//    wfGradUp[1] = 0;
//    wfGradDown[0] = 0;
//    wfGradDown[1] = 0;
    for (int i = 0; i < nParticles/2; i++)
    {
        if (k%2==0)
        {
            states[2*i]->wfGradient(grad,r[k],alpha,omega); // Should be regular indices
            grad[0] *= DSpinDownInverse[i][k/2];
            grad[1] *= DSpinDownInverse[i][k/2];
        }
        else
        {
            states[2*i+1]->wfGradient(grad,r[k],alpha,omega);
            grad[0] *= DSpinUpInverse[i][(k-1)/2];
            grad[1] *= DSpinUpInverse[i][(k-1)/2];
        }
//        states[2*i]->wfGradient(wfGradDown,r[k],alpha,omega); // Should be regular indices
//        states[2*i+1]->wfGradient(wfGradUp,r[k],alpha,omega); // CORRECT INDICES??????

//        wfGradDown[0] *= DSpinDownInverse[i][k/2];
//        wfGradDown[1] *= DSpinDownInverse[i][k/2];
//        wfGradUp[0] *= DSpinUpInverse[i][(k-1)/2];
//        wfGradUp[1] *= DSpinUpInverse[i][(k-1)/2];

//        grad[0] += wfGradUp[0] + wfGradDown[0];
//        grad[1] += wfGradUp[1] + wfGradDown[1];

//        grad[0] += wfGradUp[0]*DSpinUpInverse[i][(k-1)/2] + wfGradDown[0]*DSpinDownInverse[i][k/2]; // Indices must map to correct Spin matrices
//        grad[1] += wfGradUp[1]*DSpinUpInverse[i][(k-1)/2] + wfGradDown[1]*DSpinDownInverse[i][k/2];
//        grad[0] += wfGradUp[0]*DSpinUpInverse(i,k) + wfGradDown[0]*DSpinDownInverse(i,k);
//        grad[1] += wfGradUp[1]*DSpinUpInverse(i,k) + wfGradDown[1]*DSpinDownInverse(i,k);
    }
//    delete [] wfGradUp;
//    delete [] wfGradDown;
}

double NElectron::laplacianSlater(double **r, int k)
{
    double lap = 0;
    for (int i = 0; i < nParticles/2; i++) // Spin down
    {
        if (k%2==0)
        {
//            if(states[k]->getSpin() != -1) exit(1);
//            cout << "Spin = -1" << "   " << states[k]->getSpin() << endl;
//            lap += states[2*i]->wfLaplacian(r[k], alpha, omega) * DSpinDownInverse(i,k/2);
            lap += states[2*i]->wfLaplacian(r[k], alpha, omega) * DSpinDownInverse[i][k/2];
        }
        else
        {
//            cout << "Spin = 1" << "   " << states[k]->getSpin() << endl;;
//            lap += states[2*i+1]->wfLaplacian(r[k], alpha, omega) * DSpinUpInverse(i,(k-1)/2);
            lap += states[2*i+1]->wfLaplacian(r[k], alpha, omega) * DSpinUpInverse[i][(k-1)/2];
        }
    }
//    printf("Slater Laplacian = %10f \n", lap);
    return lap;
}

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
            psiSum += a(i,j) / (1.0/r_ij(r[i],r[j]) + beta); // a returns either a=1 or a=1/3
        }
    }
    return exp(psiSum);
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
    double commonFactor = 0;
    double r_dist = 0;
    double r_ijBeta = 0;
//    double *gradSum1 = new double[2];
//    double *gradSum2 = new double[2];
//    gradSum1[0] = 0;
//    gradSum1[1] = 0;
//    gradSum2[0] = 0;
//    gradSum2[1] = 0;

    for (int i = 0; i < nParticles; i++)
    {
//        if (i==k) continue;
//        r_dist = r_ij(r[k],r[i]);
//        commonFactor = a(k,i)/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
//        grad[0] += (r[k][0] - r[i][0]) * commonFactor;
//        grad[1] += (r[k][1] - r[i][1]) * commonFactor;
//        grad[0] += (r[i][0] - r[k][0]) * commonFactor;
//        grad[1] += (r[i][1] - r[k][1]) * commonFactor;
        if (i==k) continue;
        r_dist          = r_ij(r[i],r[k]);
        r_ijBeta        = (1.0 + beta*r_dist);
        commonFactor    = a(i,k)/(r_dist*r_ijBeta*r_ijBeta);
        grad[0]         += (r[k][0] - r[i][0])*commonFactor; // Is the difference r_i and r_k correct?
        grad[1]         += (r[k][1] - r[i][1])*commonFactor;
    }

//    printf("( %10f , %10f ) \n",grad[0],grad[1]);

//    for (int i = 0; i < k; i++) // STOP AT k-1 OR k?
//    {
//        r_dist = r_ij(r[i],r[k]);
//        commonFactor = a(i,k)/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
////        grad[0] += (r[i][0]-r[k][0])*commonFactor; // Is the difference r_i and r_k correct?
////        grad[1] += (r[i][1]-r[k][1])*commonFactor;
//        // TEMP
//        gradSum1[0] += (r[i][0]-r[k][0])*commonFactor; // i-k makes methods equal??
//        gradSum1[1] += (r[i][1]-r[k][1])*commonFactor;
//    }
//    for (int i = k+1; i < nParticles; i++)
//    {
//        r_dist = r_ij(r[k],r[i]);
//        commonFactor = a(k,i)/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
////        grad[0] += (r[k][0]-r[i][0])*commonFactor; // Is the difference r_i and r_k correct?
////        grad[1] += (r[k][1]-r[i][1])*commonFactor;
////        grad[0] += (r[i][0]-r[k][0])*commonFactor; // Is the difference r_i and r_k correct?
////        grad[1] += (r[i][1]-r[k][1])*commonFactor;
//        // TEMP
//        gradSum1[0] += (r[i][0]-r[k][0])*commonFactor;
//        gradSum1[1] += (r[i][1]-r[k][1])*commonFactor;
//    }
//    cout << "GradSum1 = ( " << std::setw(10) << gradSum1[0] << " , " << std::setw(10) << gradSum1[1] << " ) "
//         << "GradSum2 = ( " << std::setw(10) << gradSum2[0] << " , " << std::setw(10) << gradSum2[1] << " ) " << endl;
}

double NElectron::laplacianJastrow(double **r, int k, double *gradJastrow)
{
    /*
     * Returns the laplacian of the Jastrow factor.
     * Arguments:
     *  r   : particle positions
     *  k   : index of particle being moved
     */
    double lap = 0;
    double r_dist = 0;
    double commonFactor = 0;
    double r_ijBeta = 0;
    // Method 1
//    double sum1x = 0;
//    double sum1y = 0;
//    double sum2x = 0;
//    double sum2y = 0;
//    double sum3 = 0;
//    double sum4 = 0;
//    for (int i = 0; i < k; i++)  // STOP AT k-1 OR k?
//    {
//        r_dist = r_ij(r[i],r[k]);
//        r_ijBeta = (1+beta*r_dist);
//        commonFactor = a(i,k)/(r_dist*r_ijBeta*r_ijBeta);
//        sum1x += (r[i][0]-r[k][0])*commonFactor; // Is the difference r_i and r_k correct?
//        sum1y += (r[i][1]-r[k][1])*commonFactor;
//    }
//    for (int i = k+1; i < nParticles; i++)
//    {
//        r_dist = r_ij(r[k],r[i]);
//        r_ijBeta = (1+beta*r_dist);
//        commonFactor = a(k,i)/(r_dist*r_ijBeta*r_ijBeta);
//        sum2x += (r[k][0]-r[i][0])*commonFactor; // Is the difference r_i and r_k correct?
//        sum2y += (r[k][1]-r[i][1])*commonFactor;
//    }
//    lap += (sum1x-sum2x)*(sum1x-sum2x) + (sum1y-sum2y)*(sum1y-sum2y);
//    for (int i = 0; i < k; i++)
//    {
//        r_dist = r_ij(r[i],r[k]);
//        r_ijBeta = (1+beta*r_dist);
//        commonFactor = a(i,k)/(r_dist*r_ijBeta*r_ijBeta*r_ijBeta);
//        lap += (1.0-beta*r_dist)*commonFactor;
//        sum3 += (1.0-beta*r_dist)*commonFactor;
//    }
//    for (int i = k+1; i < nParticles; i++)  // STOP AT k-1 OR k?
//    {
//        r_dist = r_ij(r[k],r[i]);
//        r_ijBeta = (1+beta*r_dist);
//        commonFactor = a(k,i)/(r_dist*r_ijBeta*r_ijBeta*r_ijBeta);
//        lap += (1.0-beta*r_dist)*commonFactor;
//        sum4 += (1.0-beta*r_dist)*commonFactor;
//    }

////    printf("sum1 = ( %10f , %10f ) ", sum1x, sum1y);
////    printf("sum2 = ( %10f , %10f ) ", sum2x, sum2y);
////    printf("sum3 = %10f ", sum3);
////    printf("sum4 = %10f ", sum4);
//    printf("Lap1 = %10f ", lap);

    // Method 2
//    double sum_x = 0;
//    double sum_y = 0;
//    for (int i = 0; i < nParticles; i++)
//    {
//        if (i==k) continue;
//        r_dist = r_ij(r[i],r[k]);
//        r_ijBeta = (1.0 + beta*r_dist);
//        commonFactor = a(i,k)/(r_dist*r_ijBeta*r_ijBeta);
//        sum_x += (r[k][0] - r[i][0])*commonFactor; // Is the difference r_i and r_k correct?
//        sum_y += (r[k][1] - r[i][1])*commonFactor;
//    }
//    lap += sum_x*sum_x + sum_y*sum_y;

    lap += gradJastrow[0]*gradJastrow[0] + gradJastrow[1]*gradJastrow[1];
    for (int i = 0; i < nParticles; i++)
    {
        if (i==k) continue;
        r_dist          = r_ij(r[i],r[k]);
        r_ijBeta        = (1 + beta*r_dist);
        commonFactor    = a(i,k)/(r_dist*r_ijBeta*r_ijBeta*r_ijBeta);
        lap             += (1.0 - beta*r_dist)*commonFactor; // dim 2
//        lap += 2*commonFactor; // dim 3
    }
//    printf("Lap2 = %10f ", lap);

//    // Method 3
//    lap=0;
//    double commonFactor_i = 0;
//    double commonFactor_j = 0;
//    for (int i = 0; i < nParticles; i++)
//    {
//        if (i==k) continue;
//        r_dist = r_ij(r[k],r[i]);
//        commonFactor_i = a(k,i)/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
//        for (int j = 0; j < nParticles; j++)
//        {
//            if (j==k) continue;
//            r_dist = r_ij(r[k],r[j]);
//            commonFactor_j = a(k,j)/(r_dist*(1+beta*r_dist)*(1+beta*r_dist));
//            lap += ((r[k][0] - r[i][0])*(r[k][0] - r[j][0]) + (r[k][1] - r[i][1])*(r[k][1] - r[j][1]))*commonFactor_i*commonFactor_j;
//        }
//    }
//    for (int i = 0; i < nParticles; i++)
//    {
//        if (i==k) continue;
//        r_dist = r_ij(r[k],r[i]);
//        commonFactor_i = 1 + beta*r_dist;
//        lap += (1.0-beta*r_dist)*a(k,i)/(r_dist*commonFactor_i*commonFactor_i*commonFactor_i); // Dimension = 2
////        lap += 2*a(k,i)/(r_dist*commonFactor_i*commonFactor_i*commonFactor_i); // If dimension is d=3
//    }
//    printf("Lap3 = %10f\n", lap);
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
    double lap = 0;
    double *gradSlater = new double[2];
    double *gradJastrow = new double[2];
    gradSlater[0] = 0;
    gradSlater[1] = 0;
    gradJastrow[0] = 0;
    gradJastrow[1] = 0;
    // Running with/without Jastrow
    if (runJastrow)
    {
        gradientSlater(gradSlater,r,k);
        gradientJastrow(gradJastrow,r,k);
        lap = laplacianSlater(r,k) + laplacianJastrow(r, k, gradJastrow) + 2*(gradJastrow[0]*gradSlater[0] + gradJastrow[1]*gradSlater[1]);
    }
    else
    {
        lap = laplacianSlater(r,k);
    }

    delete [] gradSlater;
    delete [] gradJastrow;

    return lap;
}

double NElectron::a(int i, int j)
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
    if (states[i]->getSpin() + states[j]->getSpin() == 0) // Returns 1 if spins are anti-parallel
    {
        return 1.0;
    }
    else // Returns 1/3 if spins are parallel
    {
        return 0.333333333333333;
    }
}

void NElectron::revert(double **r)
{
    /*
     * Function for reverting the slater matrices-
     */
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            DSpinDown[i][j]          = DSpinDownOld[i][j];
            DSpinUp[i][j]            = DSpinUpOld[i][j];
            DSpinDownInverse[i][j]   = DSpinDownInverseOld[i][j];
            DSpinUpInverse[i][j]     = DSpinUpInverseOld[i][j];
        }
    }

}
