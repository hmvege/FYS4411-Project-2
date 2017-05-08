#include "metropolisratio.h"

#include <iostream>
MetropolisRatio::MetropolisRatio()
{
    /*
     * Default Metropolis ratio. Selects next stop from the uniform distribution
     */
}

MetropolisRatio::~MetropolisRatio()
{

}

double MetropolisRatio::Ratio(double ** rPosNew, double ** rPosOld, double newWF, double oldWF)
{
    /*
     * Default is without any importance sampling
     */
    return (newWF*newWF)/(oldWF*oldWF);
}

//double MetropolisRatio::nextStep(double oldPosition)
//{

//}
