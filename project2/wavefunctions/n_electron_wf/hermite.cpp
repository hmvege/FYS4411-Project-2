#include "hermite.h"
#include <iostream>

using std::cout;
using std::endl;

typedef double (*polyArr) (double x);

Hermite::Hermite()
{
    /*
     * Class for hermite polynomial. Hardcoded the first 11 hermite polynomials.
     */
    polynomials[0] = &H0;
    polynomials[1] = &H1;
    polynomials[2] = &H2;
    polynomials[3] = &H3;
    polynomials[4] = &H4;
    // First derivatives
    polynomialsFirstDerivative[0] = &dH0;
    polynomialsFirstDerivative[1] = &dH1;
    polynomialsFirstDerivative[2] = &dH2;
    polynomialsFirstDerivative[3] = &dH3;
    polynomialsFirstDerivative[4] = &dH4;
    // Second derivatives
    polynomialsSecondDerivative[0] = &ddH0;
    polynomialsSecondDerivative[1] = &ddH1;
    polynomialsSecondDerivative[2] = &ddH2;
    polynomialsSecondDerivative[3] = &ddH3;
    polynomialsSecondDerivative[4] = &ddH4;
}

double Hermite::getRecursive(double x, int n)
{
    /*
     * Recursion relation. Elegant, but slow.
     * Arguments:
     *  x   : function input
     *  i   : hermite rank
     */
    if (n <= 0)
    {
        return 1;
    }
    else
    {
        return 2*(x*getRecursive(x, n-1) - (((double) n)-1)*getRecursive(x, n-2));
    }
}

double Hermite::get(int n, double x)
{
    /*
     * Returns Hermite polynomial of degree n.
     */
    return polynomials[n](x);
}

double Hermite::derivative(int n, double x)
{
    /*
     * Returns the derivative of the Hermite polynomial, as given by the recursion relation:
     *  d/dx H_n(x) = 2*n*H_(n-1)(x)
     */
    return polynomialsFirstDerivative[n](x);
}

double Hermite::doubleDerivative(int n, double x)
{
    /*
     * Returns the second derivative of the Hermite polynomial, as given by the recursion relation:
     *  d^2/dx^2 H_n(x) = 4*n^2*H_(n-2)(x)
     */
    return polynomialsSecondDerivative[n](x);
}

double Hermite::getAlphaDerivative(int n, double x, double alpha, double omega)
{
    if (n == 0)
    {
        return 0.0;
    }
    else if (n == 1)
    {
        return x*sqrt(omega/alpha);
    }
    else if (n==2)
    {
        return 4*omega*x*x;
    }
    else if (n==3)
    {
        return 12*sqrt(omega*alpha)*omega*x*x*x - 6*x*sqrt(omega/alpha);
    }
}
