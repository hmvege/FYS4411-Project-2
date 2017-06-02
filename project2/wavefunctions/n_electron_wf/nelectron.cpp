#include "nelectron.h"
#include <iostream>
#include "../wavefunctions.h"
#include "hermite.h"
#include "state.h"
#include "functions.h"
#include <armadillo>
#include <mpi.h>
#include <iomanip>

using std::cout;
using std::endl;

NElectron::NElectron(int new_nParticles, int new_nDimensions, int new_numprocs, int new_processRank, double new_omega, double new_alpha, double new_beta) : WaveFunctions(new_nParticles, new_nDimensions, new_numprocs, new_processRank)
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
     * [x] implement wf functions
     * [x] implement quantum force
     * [x] implement steepest descent
     * [x] get correct results for 2 electrons->identify and remove bugs!
     * [x] fix bug that makes energy slightly too high
     * [ ] implement inverse update
     * [ ] implement ratio evaluation and storage
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
    else if (nParticles==20) maxShell = 4;
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
    // Creating a-value matrix
    a = new double*[nParticles];
    for (int i = 0; i < nParticles; i++)
    {
        a[i] = new double[nParticles];
        for (int j = 0; j < nParticles; j++)
        {
            a[i][j] = get_a(i,j);
//            cout << "a = " <<a[i][j] << " ";
//            cout << "i = " << i << endl;
//            states[i]->print();
//            cout << "j = " << j << endl;
//            states[j]->print();
//            cout << endl;
        }
    }

    // Initializing spin matrices
    DSpinDown           = arma::mat(nParticles/2,nParticles/2);
    DSpinUp             = arma::mat(nParticles/2,nParticles/2);
    DSpinDownInverse    = arma::mat(nParticles/2,nParticles/2);
    DSpinUpInverse      = arma::mat(nParticles/2,nParticles/2);
    DSpinDownOld        = arma::mat(nParticles/2,nParticles/2);
    DSpinUpOld          = arma::mat(nParticles/2,nParticles/2);
    DSpinDownInverseOld = arma::mat(nParticles/2,nParticles/2);
    DSpinUpInverseOld   = arma::mat(nParticles/2,nParticles/2);
    DSpinDown.zeros();
    DSpinDownOld.zeros();
    DSpinUp.zeros();
    DSpinUpOld.zeros();
    DSpinDownInverse.zeros();
    DSpinDownInverseOld.zeros();
    DSpinUpInverse.zeros();
    DSpinUpInverseOld.zeros();

//    // Testing
//    for (int i = 0; i < nParticles; i++)
//    {
//        cout << "i=" << i << ", State: ";
//        states[i]->print();
//    }
}

NElectron::~NElectron()
{
    for (int i = 0; i < nParticles; i++)
    {
        delete [] a[i];
    }
    delete [] a;
}

void NElectron::initializeWFSampling(double **r)
{
    /*
     * Initializing the wave function sampling.
     * Arguments:
     *  r   : position
     */
    initializeSlater(r);
}

double NElectron::initializeWaveFunction(double **r)
{
    /*
     * Returns the first setup of the wavefunction.
     * Arguments:
     *  r   : position
     */
    WFSlater = psiSlater();
    WFSlaterOld = WFSlater;
    if (runJastrow)
    {
        WFJastrow = psiJastrow(r);
        WFJastrowOld = WFJastrow;
        return WFJastrow*WFSlater;
    }
    else
    {
        return WFSlater;
    }
}

double NElectron::calculate(double **rNew, int k)
{
    /*
     * Returns the wave function for N electrons.
     * Arguments:
     *  r   : position
     *  k   : particle being moved
     */
    updateSlater(rNew, k);
    if (runJastrow)
    {
        WFJastrow = psiJastrow(rNew);
        return WFJastrow*WFSlater;
    }
    else
    {
        return WFSlater;
    }
}


void NElectron::localEnergy(double **r, double &ETotal, double &EKinetic, double &EPotential)
{
    /*
     * Returns the local energy for N electrons.
     * Arguments:
     *  r   : particle positions
     */
    double kineticEnergy = 0;
    double potentialEnergy = 0;
    for (int i = 0; i < nParticles; i++)
    {
        kineticEnergy += -0.5*laplacian(r,i);
        potentialEnergy += 0.5*omega*omega*(r[i][0]*r[i][0] + r[i][1]*r[i][1]);
    }
    if (coulombInteraction)
    {
        potentialEnergy += coulomb(r);
    }
    EKinetic = kineticEnergy;
    EPotential = potentialEnergy;
    ETotal = kineticEnergy + potentialEnergy;
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
    double *gradSlater = new double[nDimensions];
    double *gradJastrow = new double[nDimensions];
    for (int j = 0; j < nDimensions; j++)
    {
        gradSlater[j] = 0;
        gradJastrow[j] = 0;
    }
    gradientSlater(gradSlater,r,k);
    if (runJastrow) // If set to run with Jastrow
    {
        gradientJastrow(gradJastrow,r,k);
        for (int j = 0; j < nDimensions; j++)
        {
            F[k][j] = 2*(gradSlater[j] + gradJastrow[j]);
        }
    }
    else
    {
        for (int j = 0; j < nDimensions; j++)
        {
            F[k][j] = 2*gradSlater[j];
        }
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
    alpha -= SDStepLength * 2*(dPsiEAlphaSum - dPsiAlphaSum*ESum);; // Updating alpha and beta
    if (runJastrow)
    {
        beta -= SDStepLength * 2*(dPsiEBetaSum - dPsiBetaSum*ESum);
    }
    dPsiEAlphaSum = 0; // Shouldn't these be reset?
    dPsiAlphaSum = 0;
    dPsiBetaSum = 0;
    dPsiEBetaSum = 0;
}

void NElectron::sampleSD(double **r, double &E)
{
    /*
     * Sampling used by the steepest descent algorithm.
     * Arguments:
     *  r   : particle positions
     *  E   : local energy of current positions
     */
    dPsiAlpha       = alphaDerivative(r); // Derivative of WF w.r.t. alpha
    dPsiAlphaSum    += dPsiAlpha;
    dPsiEAlphaSum   += dPsiAlpha*E;
    if (runJastrow)
    {
        dPsiBeta        = betaDerivative(r); // Derivative of WF w.r.t. beta
        dPsiBetaSum     += dPsiBeta;
        dPsiEBetaSum    += dPsiBeta*E;
    }
}

void NElectron::SDStatistics(int NCycles)
{
    /*
     * Function for retrieving Steepest Descent statistics. Arguments:
     * r        : positions
     * NCycles  : Monte Carlo cycles
     */
    dPsiAlphaSum    /= double(NCycles);
    dPsiEAlphaSum   /= double(NCycles);
    if (runJastrow)
    {
        dPsiBetaSum     /= double(NCycles);
        dPsiEBetaSum    /= double(NCycles);
    }
}

void NElectron::finalizeSD()
{
    double temp_dPsiAlphaSum = 0;
    double temp_dPsiEAlphaSum = 0;
    MPI_Reduce(&dPsiAlphaSum, &temp_dPsiAlphaSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dPsiEAlphaSum, &temp_dPsiEAlphaSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    dPsiAlphaSum = temp_dPsiAlphaSum/double(numprocs);
    dPsiEAlphaSum = temp_dPsiEAlphaSum/double(numprocs);
    MPI_Bcast(&dPsiAlphaSum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dPsiEAlphaSum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (runJastrow)
    {
        double temp_dPsiBetaSum = 0;
        double temp_dPsiEBetaSum = 0;
        MPI_Reduce(&dPsiBetaSum, &temp_dPsiBetaSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&dPsiEBetaSum, &temp_dPsiEBetaSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        dPsiBetaSum = temp_dPsiBetaSum/double(numprocs);
        dPsiEBetaSum = temp_dPsiEBetaSum/double(numprocs);
        MPI_Bcast(&dPsiBetaSum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dPsiEBetaSum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

double NElectron::alphaDerivative(double **r)
{
    /*
     * The derivative of alpha w.r.t. the Slater determinant as used by steepest descent.
     * Arguments:
     * r    : particle positions
     */
    double dAlphaSpinDown = 0;
    double dAlphaSpinUp = 0;
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            dAlphaSpinDown = states[2*j]->wfAlpha(r[2*i],alpha,omega) * DSpinDownInverse(j,i);
            dAlphaSpinUp = states[2*j+1]->wfAlpha(r[2*i+1],alpha,omega) * DSpinUpInverse(j,i);
        }
    }
    return dAlphaSpinUp + dAlphaSpinDown;
}

double NElectron::betaDerivative(double **r)
{
    /*
     * The derivative of beta w.r.t. the Jastrow factor as used by steepest descent.
     * Arguments:
     * r    : particle positions
     */
    double dBeta = 0;
    double r_dist = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (i==j) continue;
            r_dist = r_ij(r[i],r[j]);
            dBeta -= a[i][j]/((1/r_dist + beta)*(1/r_dist + beta));
        }
    }
    return dBeta;
}

void NElectron::printVariationalParameters(int i)
{
    /*
     * Temporary function for printing the variational parameters used.
     */
    cout << "\nSD iteration = " << std::setw(5) << i << " Alpha = " << std::setw(10) << alpha;
    if (runJastrow)
    {
        cout << " Beta = " << std::setw(10) << beta;
    }
    cout << endl << endl;;
}

void NElectron::printUpdatedVariationalParameters()
{
    /*
     * Prints updated variational parameters from steepest descent.
     */
    cout << "Updated variational parameters: " << endl;
    cout << "Alpha = " << std::setw(10) << alpha << endl;;
    if (runJastrow)
    {
        cout << "Beta = " << std::setw(10) << beta << endl;
    }
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
     * Arguments:
     * r    : particle positions
     */
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            DSpinDown(i,j)             = states[2*j]->wf(r[2*i], alpha, omega);
            DSpinUp(i,j)               = states[2*j+1]->wf(r[2*i+1], alpha, omega);
        }
    }
//    WFSlaterOld = psiSlater();;
//    WFSlater = WFSlaterOld;
    DSpinDownInverse = arma::inv(DSpinDown);
    DSpinUpInverse = arma::inv(DSpinUp);
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            DSpinDownOld(i,j)          = DSpinDown(i,j);
            DSpinUpOld(i,j)            = DSpinUp(i,j);
            DSpinDownInverseOld(i,j)   = DSpinDownInverse(i,j);
            DSpinUpInverseOld(i,j)     = DSpinUpInverse(i,j);
        }
    }
//    if (processRank==0) {
//        printDiagnostics(r, 0);
//        exit(1);
//    }
}

void NElectron::updateSlater(double **rNew, int k)
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
            DSpinDown(i,j)  = states[2*j]->wf(rNew[2*i], alpha, omega);
            DSpinUp(i,j)    = states[2*j+1]->wf(rNew[2*i+1], alpha, omega);
        }
    }
//    if (k%2==0) {
//        for (int j = 0; j < nParticles/2; j++) // States
//        {
//            DSpinDown(k/2,j)  = states[2*j]->wf(r[k/2], alpha, omega);
//            WFSlaterOld = WFSlater;
//            WFSlater = psiSlater();
//            DSpinDownInverse = arma::inv(DSpinDown);
//        }
//    } else {
//        for (int j = 0; j < nParticles/2; j++) // States
//        {
//            DSpinUp((k-1)/2,j)    = states[2*j+1]->wf(r[(k-1)/2], alpha, omega);
//            WFSlaterOld = WFSlater;
//            WFSlater = psiSlater();
//            DSpinUpInverse = arma::inv(DSpinUp);
//        }
//    }
    DSpinDownInverse = arma::inv(DSpinDown);
    DSpinUpInverse = arma::inv(DSpinUp);
    WFSlater = psiSlater();

    if (((arma::accu(DSpinDownInverse*DSpinDown)/3.0 - 1) > 1e-10) || ((arma::accu(DSpinUpInverse*DSpinUp)/3.0 - 1) > 1e-10))
    {
        cout << "ACCU ERROR: " << arma::accu(DSpinDownInverse*DSpinDown)/3.0 << endl;
        printDiagnostics(rNew, k);
        exit(1);
    }

    if ((fabs(det(DSpinDown)) < 1e-16) || (fabs(det(DSpinUp)) < 1e-16))
    {
        printDiagnostics(rNew, k);
        exit(1);
    }
//    WFSlaterOld = WFSlater;
//    arma::inv(DSpinDownInverse,DSpinDown);
//    arma::inv(DSpinUpInverse,DSpinUp);


//    for (int i = 0; i < nParticles/2; i++)
//    {
//        for (int j = 0; j < nParticles/2; j++)
//        {
//            if (k%2==0) // Spin down
//            {
//                DSpinDownInverse(i,j) = updateInverseSlaterElement(DSpinDown, DSpinDownOld, DSpinDownInverseOld, i, j, k/2);
//            }
//            else // Spin up
//            {
//                DSpinUpInverse(i,j) = updateInverseSlaterElement(DSpinUp, DSpinUpOld, DSpinUpInverseOld, i, j, (k-1)/2);
//            }
//        }
//    }

//    double eps = 10e-15;
//    if ((fabs(arma::accu(DSpinDown*DSpinDownInverse)/(nParticles/2) - 1) > eps) || (fabs(arma::accu(DSpinUp*DSpinUpInverse)/(nParticles/2) -1) > eps)) {
//        cout << std::setprecision(18) << (arma::accu(DSpinDown*DSpinDownInverse)/(nParticles/2) - 1) << endl;
//        cout << std::setprecision(18) << (arma::accu(DSpinUp*DSpinUpInverse)/(nParticles/2) -1) << endl;
//        cout << std::setprecision(18) << DSpinDown(0,0)*DSpinDownInverse(0,0) << endl;
//        cout << std::setprecision(18) << DSpinUp(0,0)*DSpinUpInverse(0,0) << endl;
//        exit(1);
//    }

}

double NElectron::updateInverseSlaterElement(arma::mat DNew,
                                           arma::mat DOld,
                                           arma::mat DInverseOld,
                                           int k, int j, int i)
{
    /*
     * An efficient way of updating the inverse of the Slater matrices.
     * Arguments:
     *  r   : Particle positions
     *  i   : Position i in SD matrix
     *  j   : Position j in SD matrix
     *  k   : Particle that is being updated
     */
    double R = WFSlater/WFSlaterOld;
    double sum = 0;
    if (j!=i)
    {
        for (int l = 0; l < nParticles/2; l++)
        {
            sum += DNew(i,l)*DInverseOld(l,j);
        }
        return DInverseOld(k,j) - DInverseOld(k,i)/R * sum;
    }
    else
    {
        for (int l = 0; l < nParticles/2; l++)
        {
            sum += DOld(i,l)*DInverseOld(l,j);
        }
        return DInverseOld(k,i)/R * sum;
    }
}

double NElectron::psiSlater()
{
    /*
     * Returns the Slater determninant for the electrons.
     * Arguments:
     *  r   : position
     */
    return det(DSpinUp)*det(DSpinDown); // Since we divide wavefunctions on each other, we do not need factorial
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
    double *wfGrad = new double[nDimensions];
    for (int j = 0; j < nDimensions; j++)
    {
        wfGrad[j] = 0;
    }
    for (int i = 0; i < nParticles/2; i++)
    {
        if (k%2==0)
        {
            states[2*i]->wfGradient(wfGrad,r[k],alpha,omega);
            for (int j = 0; j < nDimensions; j++)
            {
                grad[j] += wfGrad[j]*DSpinDownInverse(i,k/2);
            }
        }
        else
        {
            states[2*i+1]->wfGradient(wfGrad,r[k],alpha,omega);
            for (int j = 0; j < nDimensions; j++)
            {
                grad[j] += wfGrad[j]*DSpinUpInverse(i,(k-1)/2);
            }
        }
        for (int j = 0; j < nDimensions; j++)
        {
            wfGrad[j] = 0;
        }
    }
    delete [] wfGrad;
}

double NElectron::laplacianSlater(double **r, int k)
{
    /*
     * Laplacian of the Slater determinant for a single particle.
     * Arguments:
     *  r   : particle positions
     *  k   : particle to find Laplacian of
     */
    double lap = 0;
    for (int i = 0; i < nParticles/2; i++) // Spin down
    {
        if (k%2==0)
        {
            lap += states[2*i]->wfLaplacian(r[k], alpha, omega) * DSpinDownInverse(i,k/2);
        }
        else
        {
            lap += states[2*i+1]->wfLaplacian(r[k], alpha, omega) * DSpinUpInverse(i,(k-1)/2);
        }
    }
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
            psiSum += a[i][j] / (1.0/r_ij(r[i],r[j]) + beta); // a returns either a=1 or a=1/3
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
    double commonFactor = 0;
    double r_dist = 0;
    double r_ijBeta = 0;
    for (int i = 0; i < nParticles; i++)
    {
        if (i==k) continue;
        r_dist          = r_ij(r[i],r[k]);
        r_ijBeta        = (1.0 + beta*r_dist);
        commonFactor    = a[i][k]/(r_dist*r_ijBeta*r_ijBeta);
        for (int j = 0; j < nDimensions; j++)
        {
            grad[j] += (r[k][j] - r[i][j])*commonFactor;

        }
    }
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
    lap += gradJastrow[0]*gradJastrow[0] + gradJastrow[1]*gradJastrow[1];
    for (int i = 0; i < nParticles; i++)
    {
        if (i==k) continue;
        r_dist          = r_ij(r[i],r[k]);
        r_ijBeta        = (1 + beta*r_dist);
        commonFactor    = a[i][k]/(r_dist*r_ijBeta*r_ijBeta*r_ijBeta);
        lap             += (1.0 - beta*r_dist)*commonFactor; // dim 2,
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
    double lap = 0;
    double *gradSlater = new double[nDimensions];
    double *gradJastrow = new double[nDimensions];
    for (int j = 0; j < nDimensions; j++)
    {
        gradSlater[j] = 0;
        gradJastrow[j] = 0;
    }
    if (runJastrow) // Running with/without Jastrow
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

double NElectron::get_a(int i, int j)
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
        return 0.3333333333333333;
    }
}

void NElectron::updateWF()
{
    /*
     * Function for updating internal variables in case a move is accepted.
     */
    WFSlaterOld = WFSlater;
    if (runJastrow)
    {
        WFJastrowOld = WFJastrow;
    }
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            DSpinDownOld(i,j)          = DSpinDown(i,j);
            DSpinUpOld(i,j)            = DSpinUp(i,j);
            DSpinDownInverseOld(i,j)   = DSpinDownInverse(i,j);
            DSpinUpInverseOld(i,j)     = DSpinUpInverse(i,j);
        }
    }
}

void NElectron::revert()
{
    /*
     * Function for reverting the slater matrices in case Metropolis step is rejected.
     */
    WFSlater = WFSlaterOld;
    if (runJastrow)
    {
        WFJastrow = WFJastrowOld;
    }
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            DSpinDown(i,j)          = DSpinDownOld(i,j);
            DSpinUp(i,j)            = DSpinUpOld(i,j);
            DSpinDownInverse(i,j)   = DSpinDownInverseOld(i,j);
            DSpinUpInverse(i,j)     = DSpinUpInverseOld(i,j);
        }
    }
}

void NElectron::reset()
{
    /*
     * Function for reverting the slater matrices in case Metropolis step is rejected.
     */
    WFSlater = 0;
    WFSlaterOld = 0;
    if (runJastrow)
    {
        WFJastrow = 0;
        WFJastrowOld = 0;
    }
    DSpinDown.zeros();
    DSpinDownOld.zeros();
    DSpinUp.zeros();
    DSpinUpOld.zeros();
    DSpinDownInverse.zeros();
    DSpinDownInverseOld.zeros();
    DSpinUpInverse.zeros();
    DSpinUpInverseOld.zeros();
}

std::string NElectron::getParameterString()
{
    /*
     * Returns string to be used in filename.
     */
    if (runJastrow)
    {
        return "_omega" + std::to_string(omega) + "_alpha" + std::to_string(alpha) + "_beta" + std::to_string(beta);
    }
    else
    {
        return "_omega" + std::to_string(omega) + "_alpha" + std::to_string(alpha);
    }
}

void NElectron::printDiagnostics(double **r, int k)
{
    for (int i = 0; i < 1e2; i++) cout << "=";
    cout << "\nERROR" << endl;
    cout << "Particle k" << k << endl;
    cout << "Spin down positions: " << endl;
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            cout << std::setw(10) << r[2*i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    cout << "Spin up positions: " << endl;
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            cout << std::setw(10) << r[2*i+1][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Variational parameters: " << endl;
    cout << "Alpha = " << alpha << endl;
    if (runJastrow) cout << "Beta = " << beta << endl;
    cout << "Wave functions: " << endl;
    cout << "WFSlater= " << WFSlater << endl;
    cout << "WFSlaterOld = " << WFSlaterOld << endl;
    if (runJastrow) {
        cout << "WFJastrow = " << WFJastrow << endl;
        cout << "WFJastrowOld = " << WFJastrowOld << endl;
    }
    cout << endl;
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            printf("(%2d, %2d) Spin down: %20.16f ", i, j, states[2*j]->wf(r[2*i], alpha, omega));
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < nParticles/2; i++) // Particles
    {
        for (int j = 0; j < nParticles/2; j++) // States
        {
            printf("(%2d, %2d)   Spin up: %20.16f ", i, j, states[2*j+1]->wf(r[2*i+1], alpha, omega));
        }
        cout << endl;
    }
    cout << endl;
    cout << "DSpinDown: \n" << DSpinDown    <<  endl;
    cout << "DSpinUp: \n"   << DSpinUp      <<  endl;
    cout << "DSpinDownInverse: \n" << DSpinDownInverse    <<  endl;
    cout << "DSpinUpInverse: \n"   << DSpinUpInverse      <<  endl;
    cout << endl;;
    cout << "DSpinDownOld: \n" << DSpinDownOld    <<  endl;
    cout << "DSpinUpOld: \n"   << DSpinUpOld      <<  endl;
    cout << "DSpinDownInverseOld: \n" << DSpinDownInverseOld    <<  endl;
    cout << "DSpinUpInverseOld: \n"   << DSpinUpInverseOld      <<  endl;
    cout << endl;
    for (int i = 0; i < 1e2; i++) if (processRank == 0) cout << "=";
}
