#ifndef NELECTRON_H
#define NELECTRON_H

#include "../wavefunctions.h"
#include "state.h"
#include "hermite.h"
#include <armadillo>

class NElectron : public WaveFunctions
{
private:
    // Energy
    double omega;
    // Variational parameters
    double alpha;
    double beta;
    // Wavefunction storage
    double WFJastrow;
    double WFSlater;
    double WFJastrowOld;
    double WFSlaterOld;
    // Functions used internally
    void initializeSlater(double **r);
    void updateSlater(double **r, int k);
    void updateInverseSlaterElement(double **D, double **DInverse, double **r, int i, int j, int k);
    double psiJastrow(double **r);
    double psiSlater(double **r);
    void gradientJastrow(double *grad, double **r, int k);
    void gradientSlater(double *grad, double **r, int k);
    double laplacianJastrow(double **r, int k, double *gradJastrow);
    double laplacianSlater(double **r, int k);
    double laplacian(double **r, int k);
    double a(int i, int j); // Returns value of a
    // Slater spin matrices
    double **DSpinUp;
    double **DSpinDown;
    double **DSpinUpInverse;
    double **DSpinDownInverse;
    // Storing spin-matrices in case we need to revert:
    double **DSpinUpOld;
    double **DSpinDownOld;
    double **DSpinUpInverseOld;
    double **DSpinDownInverseOld;

//    arma::mat DSpinDown;
//    arma::mat DSpinUp;
//    arma::mat DSpinDownInverse;
//    arma::mat DSpinUpInverse;
    // Matrix for the different a values (which deepends on the spin) - retrive on the fly or if-test fastest?
//    double **a;
    // Hermite class instance for use in the states(to avoid creating many different instances of the same class).
    Hermite hermite;
    // Array for storing quantum states and their wave functions
    State ** states;
    // For steepest descent
    void SDStatistics(int NCycles);
    void printVariationalParameters();
    double alphaDerivative(double **r);
    double betaDerivative(double **r);
    double dPsiAlpha = 0;
    double dPsiBeta = 0;
    double dPsiBetaSum = 0;
    double dPsiAlphaSum = 0;
    double dPsiEAlphaSum = 0;
    double dPsiEBetaSum = 0;

    // FOR RUNNING JASTROW
    bool runJastrow = true;
public:
    NElectron(int new_nParticles, int new_nDimensions, int new_nVarParams, double new_omega, double new_alpha, double new_beta);
    ~NElectron();

//    void initialize(double **r, double &WF);
    void initializeWFSampling(double **r);
    double initializeWaveFunction(double **r);
    double calculate(double **r, int k);
    double localEnergy(double **r);
    void quantumForce(double **r, double **F, int k);
    void steepestDescent(double &ESum, int NCycles);
    void sampleSD(double **r, double &E);
    bool SDConvergenceCriteria();
    void revert(double **r);
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    void setBeta(double newBeta) { beta = newBeta; }

    void setJastrow(bool jastrow) { runJastrow = jastrow; }
};

#endif // NELECTRON_H
