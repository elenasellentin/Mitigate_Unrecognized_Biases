#include "Cover.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
using namespace std;

/*
 * We derive from the base class, and implement
 * an example where a Gaussian likelihood 
 * uses the wrong covariance matrix. This is
 * an example of a biased analysis, where the
 * source of the bias is either unknown, or cannot
 * be repaired.
 */


class Example2: public Cover
{
  
private:
protected:
public:
  
    Example2(string data_vector_file, string simulations_file, int Ndata_intended, int Nsims_intended,vector<double>TrueParams, string suffix_in):Cover(data_vector_file,simulations_file,Ndata_intended, Nsims_intended,TrueParams, suffix_in)
    {
      vector<double> lowerBin(2);
      vector<double> upperBin(2);
      vector<int> spacingin(2);
      lowerBin[0] = -0.4;
      lowerBin[1] = -0.9;
      upperBin[0] = 1.8;
      upperBin[1] = 1.4;
      spacingin[0] = 100;
      spacingin[1] = 100;
      
      SetupGrid(lowerBin, upperBin, spacingin);
      AssignTheoryToGridFromFunction(-2,-2,-2,-2);
      
      TheoryAtTrueParamPoint.resize(Ndat);
      for(int i = 0; i < Ndat; i++)
      {
          if(i < 50){ TheoryAtTrueParamPoint[i] = 1.0; }
          else{  TheoryAtTrueParamPoint[i] = 0.0;   }
      }
      
      
    }
  
  
  //We still assign the correct mean
    void AssignTheoryToGridFromFunction(int p1, int p2, int gp_1, int gp_2)
    {
        //in our easy examples, we ignore the parameter arguments.
        for(int p = 0; p < pof; p++)
        {
          for(int b = p; b < pof; b++)
          {
            for(int i = 0; i < inc[p]; i++)
            {
              for(int j = 0; j < inc[b]; j++)
              {
                 TheoryAtGridPoints[p][b][i][j].resize(Ndat);
                
                 for(int z = 0; z < Ndat; z++)
                 {
                     if(z < 50)
                     {
                       TheoryAtGridPoints[p][b][i][j][z] = CoordToPhysical(p, i);
                     }
                    
                    else
                    {
                      TheoryAtGridPoints[p][b][i][j][z] = CoordToPhysical(b, j);
                      
                    }
                 }
                
              }
            }
          }
        }
    }
    
    
    void BiasPrecisionMatrix()
    {
      for(int i = 0; i < Ndat; i++)
      {
        for(int j = i+1; j < Ndat; j++)
        {
          double bias = 1e-2;
          gsl_matrix_set(invC,i,j,bias);
          gsl_matrix_set(invC,j,i,bias);
        }
      }
      
      //Testing whether still pos def matrix
      gsl_matrix* aux = gsl_matrix_alloc(Ndat, Ndat);
      gsl_matrix_memcpy(aux,invC);
      gsl_linalg_cholesky_decomp (aux);
      gsl_matrix_free(aux);
       
    }
  
  
  
};




int main()
{
  
  int ndat = 100; //100-dimensional data vector
  int nsims = 2000; //thousand 100 dimensional simulations
  
  vector<double> TrueParams(2);
  TrueParams[0] = 1;
  TrueParams[1] = 0; 
  
  
  Example2 E2("../input/Examples/DataVector.txt", "../input/Examples/Simulations.txt", 100, nsims, TrueParams, "Example2_1e-2");
  E2.BiasPrecisionMatrix();
  E2.ComputeCoverage();
  
  
  
  
  
  
  
}
