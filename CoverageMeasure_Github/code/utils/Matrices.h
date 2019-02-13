/*
 * Ln(a) Sellentin
 * Imperial College London
 * and
 * University of Geneva
 * 2017
 */


#ifndef __MATRICES_H_INCLUDED__
#define __MATRICES_H_INCLUDED__

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm> //can sort a c++ vector
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <math.h>
#include <vector>
#include <sys/stat.h>  //needed to create folders, if they don't exist
using namespace std;




void printxyz(double x, double y, double z);
void printvector(vector<double> v);
void printvector(vector<int> v);

double vec1Matvec2(gsl_vector* vec1, gsl_matrix* mat, gsl_vector* vec2, int dim);

//Matrix times vector = new vector
void Matvec(gsl_matrix* mat, gsl_vector* vec, int dim, gsl_vector* result);

void printmatrix(gsl_matrix* mat, int zeilen, int spalten);

void invmatrix(gsl_matrix* tobeinverted, unsigned int dim, gsl_matrix* inverted);

double vMv(vector<double> vec1,gsl_matrix* M, vector<double> vec2,int dim );

double chisqMat(gsl_vector* data, gsl_vector* mean, gsl_matrix* covmat, int dim);
double chisqInvMat(gsl_vector* data, gsl_vector* mean, gsl_matrix* INVcovmat, int dim);
double chisqCvecsInvmat(vector<double>data, vector<double>mean, gsl_matrix* INVcovmat);

double trace(gsl_matrix* m, int dim );
void mult4mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, gsl_matrix* m4, int dim, gsl_matrix* result);
void mult3mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, int dim, gsl_matrix* result);
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result);
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result, bool alldiag);

double determinant(gsl_matrix* f, unsigned int dim);
double lndeterminant(gsl_matrix* f, unsigned int dim) ;


bool iskrampf(double t,string callingfunction);
string strint(int i);
string strouble(double i);

void OutputFolder(string foldername);

#endif
