#ifndef STATE_H
#define STATE_H

#include "hermite.h"

class State
{
private:
    Hermite *hermite = nullptr;
    int n_x;
    int n_y;
    int spin;
public:
    State();
    ~State();
    double wf(double *r_i, double alpha, double omega);
    void *wfGradient(double *wfGrad, double *r_i, double alpha, double omega);
    double wfLaplacian(double *r_i, double alpha, double omega);

    // Getters, setters & printers
    int getNx() { return n_x; }
    int getNy() { return n_y; }
    int getSpin() { return spin; }
    void set(int new_n_x, int new_n_y, int new_spin);
    void setHermite(Hermite *newHermite) { hermite = newHermite; }
    void print();
};

#endif // STATE_H
