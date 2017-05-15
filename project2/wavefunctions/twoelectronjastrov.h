#ifndef TWOELECTRONJASTROV_H
#define TWOELECTRONJASTROV_H

#include "wavefunctions.h"

class twoElectronJastrov : public WaveFunctions
{
private:
    double omega;
    double a;
    double alpha;
    double beta;
    double C;
    double coulomb(double ** r);
public:
    twoElectronJastrov(int new_nParticles, int new_nDimensions, double new_omega, double new_alpha, double new_C, double new_a, double new_beta);

    double calculate(double **r);
    double localEnergy(double **positions);
    void quantumForce(double **positions, double **F, int k);
};

#endif // TWOELECTRONJASTROV_H
