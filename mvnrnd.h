//
//  mvnrnd.h
//  Filament
//
//  Created by Rui Ma on 1/7/15.
//  Copyright (c) 2015 Rui Ma. All rights reserved.
//
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
using namespace std;
void mvnrnd(double m[],const int n, const gsl_rng *r, double f[])
{
    gsl_matrix_view sigma = gsl_matrix_view_array(m, n, n);
    // !!!
    // note that after this line sigma.matrix is a pointer towards m
    // if the content of sigma.matrix alters, so does m
    double b[n];
    for(int i=0;i<n;i++)
        b[i]=gsl_ran_gaussian(r, 1.0); // b[i]~Gaussian(0,1) with mean 0 and std 1
    gsl_vector_view bb = gsl_vector_view_array(b, n);
    // cholesky decomposition sigma=L*L^T where L is a lower triangular matrix
    gsl_linalg_cholesky_decomp(&sigma.matrix);
    // L*bb ~ Norm(0,sigma);
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, &sigma.matrix, &bb.vector);
    //gsl_matrix_fprintf(stdout, &sigma.matrix, "%g");
    for(int i=0;i<n;i++)
        f[i]=gsl_vector_get(&bb.vector, i);
}

long randi(const gsl_rng *r, unsigned long n)
{
    return gsl_rng_uniform_int(r, n)+1;
}

void invMB(double m[],double b[],const int n, double f[])
{
    gsl_matrix_view M = gsl_matrix_view_array(m , n , n );
    gsl_vector_view B = gsl_vector_view_array(b , n );
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int c;
    gsl_linalg_LU_decomp(&M.matrix, p, &c);
    gsl_linalg_LU_solve(&M.matrix, p, &B.vector, x);
    for(int i=0;i<n;i++)
        f[i]=gsl_vector_get(x,i );
}