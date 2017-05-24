#ifndef NELECTRON_H
#define NELECTRON_H

#include "../wavefunctions.h"
#include "state.h"
#include "hermite.h"

class NElectron : public WaveFunctions
{
private:
    // Energy
    double omega;
    // Variational parameters
    double alpha;
    double beta;
    // Functions used internally
    void initializeSlater(double **r);
    void updateSlater(double **r, int k);
//    double determinant(double **A, int dim);
    double psiJastrow(double **r);
    double psiSlater(double **r);
    void gradientJastrow(double *grad, double **r, int k);
    void gradientSD(double *grad, double **r, int k);
    double laplacianJastrow(double **r, int k);
    double laplacianSD(double **r, int k);
    double laplacian(double **r, int k);
    double checkSpin(int i, int j); // Returns value of a
    // Slater spin matrices
    double **DSpinUp;
    double **DSpinDown;
    double **InverseDSpinUp;
    double **InverseDSpinDown;
    // Matrix for the different a values (which deepends on the spin) - retrive on the fly or if-test fastest?
    double **a;
    // Hermite class instance for use in the states(to avoid creating many different instances of the same class).
    Hermite hermite;
    // Array for storing quantum states and their wave functions
    State * states;
    // For steepest descent
    double dPsiAlpha = 0;
    double dPsiBeta = 0;
    double dPsiBetaSum = 0;
    double dPsiAlphaSum = 0;
    double dPsiEAlphaSum = 0;
    double dPsiEBetaSum = 0;
public:
    NElectron(int new_nParticles, int new_nDimensions, int new_nVarParams, double new_omega, double new_alpha, double new_beta);
    ~NElectron();

    void initialize(double **r);
    double calculate(double **r);
    double localEnergy(double **r);
    void quantumForce(double **r, double **F, int k);
    void steepestDescent(double &E, int NCycles);
    void sampleSD(double **r, double &E);
    bool SDConvergenceCriteria();
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    void setBeta(double newBeta) { beta = newBeta; }
};

#endif // NELECTRON_H