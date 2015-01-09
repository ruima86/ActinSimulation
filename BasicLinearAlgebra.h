//
//  BasicLinearAlgebra.h
//  Filament
//
//  Created by Rui Ma on 1/6/15.
//  Copyright (c) 2015 Rui Ma. All rights reserved.
//
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

void VecSumScalarProd(double lambda, double a[], double b[], double sum[])
{
    for(int i=0;i<3;i++)
        sum[i]=lambda*a[i]+b[i];
}

void VecScalarProd(double lambda, double a[],int n, double prod[])
{
    for(int i=0;i<n;i++)
        prod[i]=lambda*a[i];
}

void VecSet(double a[], double a1,double a2,double a3)
{
    a[0]=a1;a[1]=a2;a[2]=a3;
}

void VecPrint(double a[],int n)
{
    cout << setprecision(4);
    cout << setiosflags(ios::fixed);
    cout << setiosflags(ios::left);
    for(int i=0;i<n;i++)
        cout << "vector[" << i+1 << "]=" << a[i] << endl;        
}

void MatrixPrint(double a[],int m, int n)
{
    cout << setprecision(4);
    cout << setiosflags(ios::fixed);
    cout << setiosflags(ios::left);
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            cout << "matrix[" << i+1 << "," << j+1 << "]=" << a[i*n+j] << endl;
}

double dot(double a[], double b[])
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void VecCross(double a[], double b[], double c[])
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

double norm2(double a[])
{
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

void MatrixVecProd(double A[][3], double b[], double x[])
{
    for(int i=0;i<3;i++)
    {
        x[i]=0;
        for (int j=0; j<3; j++) {
            x[i]=x[i]+A[i][j]*b[j];
        }
    }
}

