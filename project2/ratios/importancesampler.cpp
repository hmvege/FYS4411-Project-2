#include "importancesampler.h"
#include <cmath>

ImportanceSampler::ImportanceSampler(double newD, double newDeltat, int newNPart, int newNDim)
{
    NPart = newNPart;
    NDim = newNDim;
    D = newD;
    deltat = newDeltat; // Should be between 0.001-0.01
    deltatD = D*deltat;
    exp_denom_factor = 1.0/(4*deltatD);
    denom_factor = 1.0/pow(4*M_PI*deltatD, 3*double(NPart)/2.);
}

double ImportanceSampler::Ratio(double ** rPosNew, double ** rPosOld, double newWF, double oldWF)
{
    return G(rPosNew, rPosOld)*newWF*newWF / (G(rPosOld,rPosNew)*oldWF*oldWF);
//    double **QMForce = WF->quantumForce(rPosNew);
////    double ySum = 0; // Old positions?
////    double xSum = 0; // New positions?
//    double expSquaredSum = 0;
//    double expSum = 0;
//    for (int i = 0; i < NDim; i++) {
//        for (int j = 0; j < NPart; j++) {
//            expSquaredSum += rPosOld[j][i] - rPosNew[j][i] - deltatD*QMForce[j][i];
//        }
//        expSum += expSquaredSum*expSquaredSum;
//        expSquaredSum = 0;
//    }

//    return denom_factor*exp( - expSum * exp_denom_factor);
}

double ImportanceSampler::G(double **rPosNew, double **rPosOld)
{
    /*
     * Greens function used in the importance sampling ratio.
     * rPosNew : y
     * rPosOld : x
     */
    double **QMForce = WF->quantumForce(rPosOld);
    double expSquaredSum = 0;
    double expSum = 0;
    for (int i = 0; i < NDim; i++) {
        expSquaredSum = 0;
        for (int j = 0; j < NPart; j++) {
            expSquaredSum += rPosNew[j][i] - rPosOld[j][i] - deltatD*QMForce[j][i];
        }
        expSum += expSquaredSum*expSquaredSum;
    }
    return denom_factor*exp( - expSum * exp_denom_factor);
}
