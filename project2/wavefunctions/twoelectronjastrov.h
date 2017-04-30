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
public:
    twoElectronJastrov(double omega_, double a_, double alpha_, double beta_, double C_);

    double calculate(double ** positions);
};

#endif // TWOELECTRONJASTROV_H
