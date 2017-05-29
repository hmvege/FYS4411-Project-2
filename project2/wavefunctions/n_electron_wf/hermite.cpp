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
    polynomials[5] = &H5;
    polynomials[6] = &H6;
    polynomials[7] = &H7;
    polynomials[8] = &H8;
    polynomials[9] = &H9;
    polynomials[10] = &H10;
}

double Hermite::getRecursive(double x, int n)
{
    /*
     * Recursion relation. Elegant, but slow.
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
    if (n < 11 && n >= 0) // Ensuring we do have an hermite polynomlal of degree 0<=n<=10.
    {
        return polynomials[n](x);
    }
    else
    {
        return 0.0;
    }
}

double Hermite::derivative(int n, double x)
{
    /*
     * Returns the derivative of the Hermite polynomial, as given by the recursion relation:
     *  d/dx H_n(x) = 2*n*H_(n-1)(x)
     */
    return 2*n*get(n-1, x);
}

double Hermite::doubleDerivative(int n, double x)
{
    /*
     * Returns the second derivative of the Hermite polynomial, as given by the recursion relation:
     *  d^2/dx^2 H_n(x) = 4*n^2*H_(n-2)(x)
     */
    return 2*n*derivative(n-1,x);
}
