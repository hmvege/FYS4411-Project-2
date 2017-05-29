#include "functions.h"

using std::cout;
using std::endl;
using std::nothrow;

int factorial(int n)
{
    /*
     * Returns the factorial of n.
     */
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double determinant(double **A, int dim)
{
    /*
     * Function for calculating the determinant.
     * Taken from lecture.
     * Arguments:
     *  A   : square matrix
     *  dim : dimensionality of the matrix
     */
    if (dim == 1) // Is this safe?
    {
        return A[0][0];
    }
    if (dim == 2)
    {
        return A[0][0]*A[1][1] - A[0][1]*A[1][0];
    }
    double sum = 0;
    for (int i = 0; i < dim; i++)
    {
        double ** sub = new double*[dim-1];
        for (int j = 0; j < i; j++)
        {
            sub[j] = &A[j][1];
        }
        for (int j = i + 1; j < dim; j++)
        {
            sub[j-1] = &A[j][1];
        }
        if (i % 2 == 0)
        {
            sum += A[i][0] * determinant(sub, dim-1);
        }
        else
        {
            sum -= A[i][0] * determinant(sub, dim-1);
        }
    }
    return sum;
}

/*
 * The function
 *      void  **matrix()
 * reserves dynamic memory for a two-dimensional matrix
 * using the C++ command new . No initialization of the elements.
 * Input data:
 *  int row      - number of  rows
 *  int col      - number of columns
 *  int num_bytes- number of bytes for each
 *                 element
 * Returns a void  **pointer to the reserved memory location.
 */

void **matrix(int row, int col, int num_bytes)
{
    int      i, num;
    char     **pointer, *ptr;

    pointer = new(nothrow) char* [row];
    if(!pointer) {
        cout << "Exception handling: Memory allocation failed";
        cout << " for "<< row << "row addresses !" << endl;
        return NULL;
    }
    i = (row * col * num_bytes)/sizeof(char);
    pointer[0] = new(nothrow) char [i];
    if(!pointer[0]) {
        cout << "Exception handling: Memory allocation failed";
        cout << " for address to " << i << " characters !" << endl;
        return NULL;
    }
    ptr = pointer[0];
    num = col * num_bytes;
    for(i = 0; i < row; i++, ptr += num )   {
        pointer[i] = ptr;
    }

    return  (void **)pointer;

} // end: function void **matrix()


/*
 * The function
 *      void free_matrix()
 * releases the memory reserved by the function matrix()
 *for the two-dimensional matrix[][]
 * Input data:
 *  void far **matr - pointer to the matrix
 */

void free_matrix(void **matr)
{

    delete [] (char *) matr[0];
    delete [] matr;

}  // End:  function free_matrix()

void inverse(double **a, int n)
{
    /*
     * Returns the matrix a as an inverse
     */
    int          i,j, *indx;
    double       d, *col, **y;

    // allocate space in memory
    indx = new int[n];
    col  = new double[n];
    y    = (double **) matrix(n, n, sizeof(double));

    ludcmp(a, n, indx, &d);   // LU decompose  a[][]

    // find inverse of a[][] by columns

    for(j = 0; j < n; j++) {

        // initialize right-side of linear equations

        for(i = 0; i < n; i++) col[i] = 0.0;
        col[j] = 1.0;

        lubksb(a, n, indx, col);

        // save result in y[][]

        for(i = 0; i < n; i++) y[i][j] = col[i];

    }   //j-loop over columns

    // return the inverse matrix in a[][]

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) a[i][j] = y[i][j];
    }
    free_matrix((void **) y);     // release local memory
    delete [] col;
    delete []indx;

}  // End: function inverse()

/*
** The function
**       ludcmp()
** takes as input a two-dimensional matrix a[][] of dimension n and
** replaces it by the LU decomposition of a rowwise permutation of
** itself. The results is stored in a[][] in the form given by
** eq. (2.3.14) in "Numerical Recipe", sect. 2.3, page 45. The vector
** indx[] records the row permutation effected by the partial pivoting;
** d is output as +1 or -1 depending on whether the number of row
** interchanges was even or odd, respectively. This routine is used in
** combination with the function lubksb() to solve linear equations or
** invert a matrix. The function is slightly modified from the version
** in in Numerical recipe and uses memory allocation functions in the
** present module.
*/

void ludcmp(double **a, int n, int *indx, double *d)
{
    int      i, imax, j, k;
    double   big, dum, sum, temp, *vv;

    vv = new(nothrow) double [n];
    if(!vv) {
        printf("\n\nError in function ludcm():");
        printf("\nNot enough memory for vv[%d]\n",n);
        exit(1);
    }

    *d = 1.0;                              // no row interchange yet
    for(i = 0; i < n; i++) {     // loop over rows to get scaling information
        big = ZERO;
        for(j = 0; j < n; j++) {
            if((temp = fabs(a[i][j])) > big) big = temp;
        }
        if(big == ZERO) {
            printf("\n\nSingular matrix in routine ludcmp()\n");
//            exit(1);
        }
        vv[i] = 1.0/big;                 // save scaling */
    } // end i-loop */

    for(j = 0; j < n; j++) {     // loop over columns of Crout's method
        for(i = 0; i< j; i++) {   // not i = j
            sum = a[i][j];
            for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
        }
        big = ZERO;   // initialization for search for largest pivot element
        for(i = j; i< n; i++) {
            sum = a[i][j];
            for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if((dum = vv[i]*fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        } // end i-loop
        if(j != imax) {    // do we need to interchange rows ?
            for(k = 0;k< n; k++) {       // yes
                dum        = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k]    = dum;
            }
            (*d)    *= -1;            // and change the parit of d
            vv[imax] = vv[j];         // also interchange scaling factor
        }
        indx[j] = imax;
        if(fabs(a[j][j]) < ZERO)  a[j][j] = ZERO;

        /*
    ** if the pivot element is zero the matrix is singular
    ** (at least to the precision of the algorithm). For
    ** some application of singular matrices, it is desirable
    ** to substitute ZERO for zero,
    */

        if(j < (n - 1)) {                   // divide by pivot element
            dum = 1.0/a[j][j];
            for(i=j+1;i < n; i++) a[i][j] *= dum;
        }
    } // end j-loop over columns

    delete [] vv;   // release local memory

}  // End: function ludcmp()

/*
** The function
**             lubksb()
** solves the set of linear equations A X = B of dimension n.
** a[][] is input, not as the matrix A[][] but rather as
** its LU decomposition, determined by the function ludcmp(),
** indx[] is input as the permutation vector returned by
** ludcmp(). b[] is input as the right-hand side vector B,
** The solution X is returned in B. The input data a[][],
** n and indx[] are not modified. This routine take into
** account the possibility that b[] will begin with many
** zero elements, so it is efficient for use in matrix
** inversion.
** The function is slightly modified from the version in
** in Numerical recipe.
*/

void lubksb(double **a, int n, int *indx, double *b)
{
    int        i, ii = -1, ip, j;
    double     sum;

    for(i = 0; i< n; i++) {
        ip    = indx[i];
        sum   = b[ip];
        b[ip] = b[i];
        if(ii > -1)   for(j = ii; j < i; j++) sum -= a[i][j] * b[j];
        else if(sum) ii = i;
        b[i] = sum;
    }
    for(i = n - 1; i >= 0; i--) {
        sum = b[i];
        for(j = i+1; j < n; j++) sum -= a[i][j] * b[j];
        b[i] = sum/a[i][i];
    }
} // End: function lubksb()

void printMatrix(double **A, int dim)
{
    /*
     * Simple matrix printer.
     */
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
//            printf("%10.5f ", A[i][j]);
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

void multiplieMatrices(double **A, double **B, double **C, int dim)
/*
 * Simple matrix multiplier.
 */
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            for (int k = 0; k < dim; k++)
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}
