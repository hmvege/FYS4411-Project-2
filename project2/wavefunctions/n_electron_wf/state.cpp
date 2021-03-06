#include "state.h"
#include "hermite.h"
#include <iostream>
#include <cmath>

State::State()
{

}

State::~State()
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
    return hermite->get(n_x, sqrtOmegaAlpha*r_i[0]) * hermite->get(n_y, sqrtOmegaAlpha*r_i[1]) * exp(-alpha*omega*0.5*(r_i[0]*r_i[0] + r_i[1]*r_i[1]));
//    return exp(-alpha*omega*0.5*(r_i[0]*r_i[0] + r_i[1]*r_i[1])); // plain and simple omega
}

void State::wfGradient(double * wfGrad, double *r_i, double alpha, double omega)
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
    double HermXDerivative = sqrtOmegaAlpha*hermite->derivative(n_x, sqrtOmegaAlpha*r_i[0]);
    double HermYDerivative = sqrtOmegaAlpha*hermite->derivative(n_y, sqrtOmegaAlpha*r_i[1]);
    double expFactor = exp(-alpha*omega*0.5*(r_i[0]*r_i[0] + r_i[1]*r_i[1]));
    wfGrad[0] = (HermXDerivative - omega*alpha*r_i[0]*HermX)*HermY*expFactor;
    wfGrad[1] = (HermYDerivative - omega*alpha*r_i[1]*HermY)*HermX*expFactor;
}

double State::wfLaplacian(double *r_i, double alpha, double omega)
{
    /*
     * Returns the gradient of the wave function. Index 0 is x, 1 is y.
     * Arguments:
     * r_i      : (x,y) position of particle i
     * alpha    : variational paramater
     * omega    : oscillator frequency
     */
    double wfLap;
    double omegaAlpha = omega*alpha;
    double sqrtOmegaAlpha = sqrt(omegaAlpha);
    double HermX = hermite->get(n_x, sqrtOmegaAlpha*r_i[0]);
    double HermY = hermite->get(n_y, sqrtOmegaAlpha*r_i[1]);
    double HermXDerivative = sqrtOmegaAlpha*hermite->derivative(n_x, sqrtOmegaAlpha*r_i[0]);
    double HermYDerivative = sqrtOmegaAlpha*hermite->derivative(n_y, sqrtOmegaAlpha*r_i[1]);
    double HermXSecondDerivative = omegaAlpha*hermite->doubleDerivative(n_x,sqrtOmegaAlpha*r_i[0]);
    double HermYSecondDerivative = omegaAlpha*hermite->doubleDerivative(n_y,sqrtOmegaAlpha*r_i[1]);
    double rr = r_i[0]*r_i[0] + r_i[1]*r_i[1];
    double expFactor = exp(-omegaAlpha*0.5*rr);
//    wfLap = ( (HermY*hermite->doubleDerivative(n_x,sqrtOmegaAlpha*r_i[0])+ HermX*hermite->doubleDerivative(n_y,sqrtOmegaAlpha*r_i[1]))
//            - 2*omegaAlpha*(r_i[0]*HermY*HermXDerivative + r_i[1]*HermX*HermYDerivative)
//            - 2*omegaAlpha*HermX*HermY + omegaAlpha*omegaAlpha*(r_i[0]*r_i[0] + r_i[1]*r_i[1])*HermX*HermY)*expFactor;
//    printf("d^2H/dx^2 = %20.16f  d^2H/dx^2 = %20.16'f \n", HermXDerivative, HermYDerivative);
    wfLap = ((HermY*HermXSecondDerivative + HermX*HermYSecondDerivative)
            - 2*omegaAlpha*(r_i[0]*HermY*HermXDerivative + r_i[1]*HermX*HermYDerivative)
            + omegaAlpha*(omegaAlpha*rr - 2)*HermX*HermY)*expFactor; // REWRITE SUCH THAT LITTLE OR NO MEMORY IS WRITTEN!
    return wfLap;
}

double State::wfAlpha(double *r_i, double alpha, double omega)
{
    /*
     * The derivative of phi w.r.t. variational parameter alpha.
     * Arguments:
     * r_i      : (x,y) position of particle i
     * alpha    : variational paramater
     * omega    : oscillator frequency
     */
    double wfA = 0;
    double sqrtOmegaAlpha = sqrt(alpha*omega);
    double HermX = hermite->get(n_x, sqrtOmegaAlpha*r_i[0]);
    double HermY = hermite->get(n_y, sqrtOmegaAlpha*r_i[1]);
    double expFactor = exp(-alpha*omega*0.5*(r_i[0]*r_i[0] + r_i[1]*r_i[1]));
    wfA = ( sqrt(omega/alpha) * (  r_i[0]*HermY*hermite->derivative(n_x, sqrtOmegaAlpha*r_i[0])
                                 + r_i[1]*HermX*hermite->derivative(n_y, sqrtOmegaAlpha*r_i[1]))
                  - omega*(r_i[0]*r_i[0] + r_i[1]*r_i[1])*HermX*HermY)*0.5*expFactor;
//    wfA = ( (  HermY*hermite->getAlphaDerivative(n_x, r_i[0], alpha, omega)
//             + HermX*hermite->getAlphaDerivative(n_y, r_i[1], alpha, omega))
//                  - 0.5*omega*(r_i[0]*r_i[0] + r_i[1]*r_i[1])*HermX*HermY)*expFactor;
    return wfA;
}
