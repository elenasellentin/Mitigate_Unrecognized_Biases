/*
 * Ln(a) Sellentin
 * Imperial College London
 * and
 * University of Geneva
 * 2017
 */

#include "Matrices.h"
using namespace std;

void printxyz(double x, double y, double z)
{
    printf("%.3f %.3f %.3f \n",x,y,z);
}


void printvector(vector<double> v)
{
    for(int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " " << endl;
    }
}


void printvector(vector<int> v)
{
    for(int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " " << endl;
    }
}


double vec1Matvec2(gsl_vector* vec1, gsl_matrix* mat, gsl_vector* vec2, int dim)
{
    gsl_vector * hilf = gsl_vector_alloc(dim);
    double scalarResult;
    gsl_blas_dgemv(CblasNoTrans, 1.0, mat, vec2, 0.0, hilf);  //Mat*vec2
    gsl_blas_ddot(hilf,vec1 , &scalarResult);

    gsl_vector_free(hilf);
    return scalarResult;
}

void Matvec(gsl_matrix* mat, gsl_vector* vec, int dim, gsl_vector* result)
{
    gsl_blas_dgemv(CblasNoTrans, 1.0, mat, vec, 0.0, result);  //Mat*vec2
}







void printmatrix(gsl_matrix* mat, int zeilen, int spalten)
{
    for(int i = 0; i < zeilen; i ++)
    {
        for(int j = 0; j < spalten; j++)
        {
            printf("%.2e ",gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
}


void invmatrix(gsl_matrix* tobeinverted, unsigned int dim, gsl_matrix* inverted)
{
    /*Takes the matrix "tobeinverted", which is of dimension "dim", and inverts it. The result is stored in an externally provided gsl_matrix "inverted". */

    /*Do not allocate space for the resulting matrix within this function!
     The code would compile, and store the matrix here - but you want to use it outside this function. However, the outside allocated matrix will then not be the one you just calculated!*/
    gsl_matrix* hilf = gsl_matrix_alloc(dim,dim);

    /*use memcopy, such that the content of the matrix is copied.
     * Else, the pointers would be set equal.*/
    gsl_matrix_memcpy(hilf,tobeinverted);
    int signum;

    gsl_permutation * perm = gsl_permutation_alloc(dim);
    gsl_linalg_LU_decomp(hilf, perm, &signum);
    gsl_linalg_LU_invert(hilf,perm,inverted);

    gsl_matrix_free(hilf);
    gsl_permutation_free(perm);
    /* cout << "Inverted Matrix:" << endl;
     printmatrix(inverted,dim,dim);*/

}






double determinant(gsl_matrix* f, unsigned int dim)
{
    /*returns the determinant of the square matrix f, when dim is the dimension of f*/

    gsl_matrix* hilf = gsl_matrix_alloc(dim,dim);
    gsl_matrix_memcpy(hilf,f); //use memcopy, to not delete the external matrix f
    int signum;

    gsl_permutation * perm = gsl_permutation_alloc(dim);
    gsl_linalg_LU_decomp(hilf, perm, &signum);
    double det = gsl_linalg_LU_det(hilf,signum);

    gsl_permutation_free(perm);
    gsl_matrix_free(hilf);
    return det;
}



double lndeterminant(gsl_matrix* f, unsigned int dim) //only works for posdef matrices
{
    /*returns the determinant of the square matrix f, when dim is the dimension of f*/

    gsl_matrix* hilf = gsl_matrix_alloc(dim,dim);
    gsl_matrix_memcpy(hilf,f); //use memcopy, to not delete the external matrix f

    /* double scale = 1./(0.25*( gsl_matrix_get(hilf,0,0 ) + gsl_matrix_get(hilf,(int)(dim/3), (int)(dim/3) )
                         + gsl_matrix_get(hilf,(int)(2*dim/3),(int)(2*dim/3) ) + gsl_matrix_get(hilf,dim-1, dim-1 )   )
                        );*/

    // gsl_matrix_scale(hilf,scale);



    int signum;

    gsl_permutation * perm = gsl_permutation_alloc(dim);
    gsl_linalg_LU_decomp(hilf, perm, &signum);

    //double lndet = log(gsl_linalg_LU_det(hilf,signum) )  - (double)(dim)*log(scale) ;

    double lndet = gsl_linalg_LU_lndet(hilf); //returns log ( fabs(det) ), therefore, no neg. def. matrices
    //double det = exp(lndet);
    gsl_permutation_free(perm);

    /* double det = 1.0;
     for(int i = 0; i < dim; i++)
     {
       det*=gsl_matrix_get(hilf,i,i);
       cout << det << endl;
     }
     cout << endl;*/


    gsl_matrix_free(hilf);
    iskrampf(lndet, __func__);
    return lndet;
}




/*scalar product without gsl_vector*/
/*assuming M is square with size dim times dim*/
double vMv(vector<double> vec1,gsl_matrix* M, vector<double> vec2,int dim )
{
    double res = 0;

    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            res+=vec1[i]*gsl_matrix_get(M, i, j)*vec2[j];
        }
    }


    return res;
}


double trace(gsl_matrix* m, int dim )
{
    double tr = 0.0;

    for(int i = 0; i < dim; i++)
    {
        tr += gsl_matrix_get(m,i,i);
    }

    return tr;
}



double chisqMat(gsl_vector* data, gsl_vector* mean, gsl_matrix* covmat, int dim)
{
    gsl_matrix* invC = gsl_matrix_calloc(dim, dim);
    invmatrix(covmat, dim, invC);
    gsl_vector* xmimu = gsl_vector_calloc(dim);

    for(int i = 0; i < dim; i++)
    {
        double temp = gsl_vector_get(data,i) - gsl_vector_get(mean, i);
        gsl_vector_set(xmimu, i, temp);
    }

    double chi2 = vec1Matvec2(xmimu, invC, xmimu, dim);

    gsl_vector_free(xmimu);
    gsl_matrix_free(invC);

    return chi2;
}



double chisqInvMat(gsl_vector* data, gsl_vector* mean, gsl_matrix* INVcovmat, int dim)
{

    gsl_vector* xmimu = gsl_vector_calloc(dim);

    for(int i = 0; i < dim; i++)
    {
        double temp = gsl_vector_get(data,i) - gsl_vector_get(mean, i);
        gsl_vector_set(xmimu, i, temp);
    }

    double chi2 = vec1Matvec2(xmimu, INVcovmat, xmimu, dim);

    gsl_vector_free(xmimu);

    return chi2;
}





double chisqCvecsInvmat(vector<double>data, vector<double>mean, gsl_matrix* INVcovmat)
{

   int dim = data.size();
   if(dim < 0.0)
   {
     cout << "Negative dimension in File Matrices.cpp, function chisqCvecsInvmat" << endl;
     cout << dim << endl;
   }
   

   
   gsl_vector* gsl_data = gsl_vector_calloc(dim);
   gsl_vector* gsl_mean = gsl_vector_calloc(dim);
   
   for(int i=0; i < dim; i++)
   {
     gsl_vector_set(gsl_data,i,data[i]);   
     gsl_vector_set(gsl_mean,i,mean[i]);   
   }
    
   double r = chisqInvMat(gsl_data, gsl_mean, INVcovmat, dim);
    
   gsl_vector_free(gsl_data);
   gsl_vector_free(gsl_mean);
   
   return r;
}







/*multiplies four matrices*/
void mult4mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, gsl_matrix* m4, int dim, gsl_matrix* result)
{

    gsl_matrix* intermed1 = gsl_matrix_calloc(dim,dim);
    gsl_matrix* intermed2 = gsl_matrix_calloc(dim,dim);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m3,m4,0.0,intermed1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m2,intermed1,0.0,intermed2);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,intermed2,0.0,result);

    gsl_matrix_free(intermed1);
    gsl_matrix_free(intermed2);
}



void mult3mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, int dim, gsl_matrix* result)
{

    gsl_matrix* intermed1 = gsl_matrix_calloc(dim,dim);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m2,m3,0.0,intermed1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,intermed1,0.0,result);

    gsl_matrix_free(intermed1);
}

/*multiplying two matrices, and storing the result in result. The matrix result must be *externally* allocated.*/
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result)
{
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,m2,0.0,result);

}

/*flag alldiag shall be true if m1 and m2 are diagonal matrices*/
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result, bool alldiag)
{
    if(!alldiag)
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,m2,0.0,result);


    if(alldiag)
    {
        gsl_matrix_set_identity(result);
        for(int i = 0; i < dim; i++)
        {
            gsl_matrix_set(result,i,i,gsl_matrix_get(m1,i,i)*gsl_matrix_get(m2,i,i));
        }
    }

}






/**returns yes, if argument is krampf;
 * if one puts __func__ into the string,
 * then it prints out the name of the
 * function that produced the nonesense
 **/
bool iskrampf(double t, string callingfunction)
{
    bool yes = false;
    if( isnan(t) ) yes = true;
    if( isinf(t) ) yes = true;

    if(yes) cout << "Isnan or isinf in function " << callingfunction << endl;

    return yes;
}


/**Takes an integer and returns it as a
 * string; handy to produce automatic
 * filenames
 **/
string strint(int i)
{
    stringstream I;
    I << i;
    return I.str();
}


string strouble(double i)
{
    stringstream I;
    I << i;
    return I.str();
}



/**Creates a new outputfolder, if it doesn't exist.**/
void OutputFolder(string foldername)
{
    /*Function returns zero, if directory exists*/
    struct stat dircheck;
    if(stat(foldername.c_str(),&dircheck) == 0 && S_ISDIR(dircheck.st_mode))
    {
        //cout << "#Output directory: " << outdir << endl;
    }

    else /*Make directory, if not existent*/
    {

        string createdirectory = (string)"mkdir " + foldername;
        int m = system(createdirectory.c_str());
        cout << "#Newly created output directory: " << foldername << endl;
    }
}

