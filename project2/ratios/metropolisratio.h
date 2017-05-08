#ifndef METROPOLISRATIO_H
#define METROPOLISRATIO_H


class MetropolisRatio
{
public:
    MetropolisRatio();
    virtual ~MetropolisRatio();

    virtual double Ratio(double ** rPosNew, double ** rPosOld, double newWF, double oldWF);
//    virtual double nextStep();
};

#endif // METROPOLISRATIO_H
