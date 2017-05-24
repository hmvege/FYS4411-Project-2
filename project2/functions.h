#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <iostream>

#define   NULL_PTR   (void *) 0
#define   ZERO       1.0E-10
#define   UL         unsigned long

int factorial(int n);
double determinant(double **A, int dim);
void inverse(double **a, int n);
void free_matrix(void **);
void **matrix(int row, int col, int num_bytes);
void ludcmp(double **, int, int *, double*);
void lubksb(double **, int, int *, double *);


#endif // FUNCTIONS_H
