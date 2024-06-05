#ifndef LAMBDA_H
#define LAMBDA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <Eigen/Dense>

using namespace Eigen;
#define MAXCOL      1024
#define LOOPMAX     100000           /* maximum count of search loop */
#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

extern double* mat(int n, int m);

int* imat(int n, int m);

double* zeros(int n, int m);

double* eye(int n);

void matmul(const char* tr, int n, int k, int m, double alpha, const double* A, const double* B, double beta, double* C);

void matcpy(double* A, const double* B, int n, int m);

int solve(const char* tr, const double* A, const double* Y, int n, int m, double* X);

static int ludcmp(double* A, int n, int* indx, double* d);

static void lubksb(const double* A, int n, const int* indx, double* b);

int matinv(double* A, int n);

static int LD(int n, const double* Q, double* L, double* D);

static int LDLT(int n, const double* Q, double* L, double* D);

static void gauss(int n, double* L, double* Z, int i, int j);

static void perm(int n, double* L, double* D, int j, double del, double* Z);

static void reduction(int n, double* L, double* D, double* Z);

static int search(int n, int m, const double* L, const double* D, const double* zs, double* zn, double* s);

/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions , m=2
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residuals of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double* a, const double* Q, double* F, double* s);

// 将MatrixXd转换为double*
double* matrixXdToDoublePtr(const MatrixXd & matrix);

// 将VectorXd转换为double*
double* vectorXdToDoublePtr(const VectorXd & vector);

// 将double*转换为MatrixXd
MatrixXd doublePtrToMatrixXd(const double* ptr, int rows, int cols);
MatrixXd doublePtrToMatrixXd1(const double* ptr, int rows, int cols);
// 将double*转换为VectorXd
VectorXd doublePtrToVectorXd(const double* ptr, int size);


#endif