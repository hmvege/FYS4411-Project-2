#include "state.h"
#include "hermite.h"
#include <iostream>
#include <cmath>

State::State()
{

}

void State::set(int new_n_x, int new_n_y, int new_spin)
{
    /*
     * For initializing parameters.
     */
    n_x = new_n_x;
    n_y = new_n_y;
    spin = new_spin;
}

void State::print()
{
    /*
     * Small function that prints the quantum numbers of the system.
     */
    printf("n_x = %3d, n_y = %3d, spin = %3d\n", n_x, n_y, spin);
}

double State::wf(double *r_i, double alpha, double omega)
{
    /*
     * Returns wavefunction at a given value.
     * Arguments
     * r_i      : (x,y) position of particle i
     * alpha    : variational paramater
     * omega    : oscillator frequency
     */
    double sqrtOmegaAlpha = sqrt(omega*alpha);
    return hermite->get(n_x, sqrtOmegaAlpha*r_i[0])*hermite->get(n_y, sqrtOmegaAlpha*r_i[1])*exp(-alpha*omega*0.5*(r_i[0]*r_i[0] + r_i[1]*r_i[1]));
}

void *State::wfGradient(double * wfGrad, double *r_i, double alpha, double omega)
{
    /*
     * Returns the gradient of the wave function. Index 0 is x, 1 is y.
     * Arguments
     * wfGrad   : gradient, pass by reference
     * r_i      : (x,y) position of particle i
     * alpha    : variational paramater
     * omega    : oscillator frequency
     */
    double sqrtOmegaAlpha = sqrt(alpha*omega);
    double HermX = hermite->get(n_x, sqrtOmegaAlpha*r_i[0]);
    double HermY = hermite->get(n_y, sqrtOmegaAlpha*r_i[1]);
    double expFactor = exp(-alpha*omega*0.5*(r_i[0]*r_i[0] + r_i[1]*r_i[1]));
    wfGrad[0] = (hermite->derivative(n_x, sqrtOmegaAlpha*r_i[0]) - omega*alpha*r_i[0]*HermX)*HermY*expFactor;
    wfGrad[1] = (hermite->derivative(n_y, sqrtOmegaAlpha*r_i[1]) - omega*alpha*r_i[1]*HermY)*HermX*expFactor;
}

double State::wfLaplacian(double *r_i, double alpha, double omega)
{
    /*
     * Returns the gradient of the wave function. Index 0 is x, 1 is y.
     * Arguments
     * r_i      : (x,y) position of particle i
     * alpha    : variational paramater
     * omega    : oscillator frequency
     */
    double wfLap;
    double omegaAlpha = omega*alpha;
    double sqrtOmegaAlpha = sqrt(omegaAlpha);
    double HermX = hermite->get(n_x, sqrtOmegaAlpha*r_i[0]);
    double HermY = hermite->get(n_y, sqrtOmegaAlpha*r_i[1]);
    double HermXDerivative = hermite->derivative(n_x, sqrtOmegaAlpha*r_i[0]);
    double HermYDerivative = hermite->derivative(n_y, sqrtOmegaAlpha*r_i[1]);
    double expFactor = exp(-omegaAlpha*0.5*(r_i[0]*r_i[0] + r_i[1]*r_i[1]));
    wfLap = ( (HermY*2*n_x*HermXDerivative + HermX*2*n_y*HermYDerivative)
            - 2*omegaAlpha*(r_i[0]*HermY*HermXDerivative + r_i[1]*HermX*HermYDerivative)
            - 2*omegaAlpha*HermX*HermY + omegaAlpha*omegaAlpha*(r_i[0]*r_i[0] + r_i[1]*r_i[1])*HermX*HermY)*expFactor;
    return wfLap;
}