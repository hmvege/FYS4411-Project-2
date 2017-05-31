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
    double updateInverseSlaterElement(arma::mat DNew, arma::mat DOld, arma::mat DInverseOld, int k, int j, int i);
    double psiJastrow(double **r);
    double psiSlater();
    void gradientJastrow(double *grad, double **r, int k);
    void gradientSlater(double *grad, double **r, int k);
    double laplacianJastrow(double **r, int k, double *gradJastrow);
    double laplacianSlater(double **r, int k);
    double laplacian(double **r, int k);
    double get_a(int i, int j); // Returns value of a
    double **a; // Matrix for the different a values(they depend on spin
    // Slater spin matrices
    arma::mat DSpinUp;
    arma::mat DSpinDown;
    arma::mat DSpinUpInverse;
    arma::mat DSpinDownInverse;
    // Storing spin-matrices in case we need to revert:
    arma::mat DSpinUpOld;
    arma::mat DSpinDownOld;
    arma::mat DSpinUpInverseOld;
    arma::mat DSpinDownInverseOld;
    // Hermite class instance for use in the states(to avoid creating many different instances of the same class).
    Hermite hermite;
    // Array for storing quantum states and their wave functions
    State ** states;
    // For steepest descent
    void SDStatistics(int NCycles);
    void printVariationalParameters(int i);
    double alphaDerivative(double **r);
    double betaDerivative(double **r);
    double dPsiAlpha = 0;
    double dPsiBeta = 0;
    double dPsiBetaSum = 0;
    double dPsiAlphaSum = 0;
    double dPsiEAlphaSum = 0;
    double dPsiEBetaSum = 0;
    bool runJastrow;
public:
    using WaveFunctions::localEnergy;
    NElectron(int new_nParticles, int new_nDimensions, int new_numprocs, int new_processRank, double new_omega, double new_alpha, double new_beta);
    ~NElectron();

//    void initialize(double **r, double &WF);
    void initializeWFSampling(double **r);
    double initializeWaveFunction(double **r);
    double calculate(double **r, int k);
    virtual void localEnergy(double **r, double &ETotal, double &EKinetic, double &EPotential);
    void quantumForce(double **r, double **F, int k);
    void steepestDescent(double &ESum, int NCycles);
    void sampleSD(double **r, double &E);
    bool SDConvergenceCriteria();
    void revert();
    void updateWF();
    std::string getParameterString();
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    void setBeta(double newBeta) { beta = newBeta; }

    void setJastrow(bool jastrow) { runJastrow = jastrow; }
};

#endif // NELECTRON_H
