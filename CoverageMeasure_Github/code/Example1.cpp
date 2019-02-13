#include "Cover.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;

/*
 * We derive from the base class, and implement
 * an example where a Gaussian likelihood 
 * times a flat prior produces the correct coverage,
 * if the means itself are inferred.
 */


class Example1: public Cover
{
  
private:
protected:
public:
  
    Example1(string data_vector_file, string simulations_file, int Ndata_intended, int Nsims_intended,vector<double>TrueParams, string suffix_in):Cover(data_vector_file,simulations_file,Ndata_intended, Nsims_intended,TrueParams, suffix_in)
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
  
  
  
};









int main()
{
  
  /*
   * Section 1: We produce 
   * - the simulations
   * - the data  for this example.
   */
  
  int ndat = 100; //100-dimensional data vector
  int nsims = 2000; //thousand 100 dimensional simulations
  
  gsl_rng * r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set (r, 11011);
  
  ofstream Datfile;
  ofstream Covfile;
  ofstream Simfile;
  
  Datfile.open("../input/Examples/DataVector.txt");
  //Covfile.open("../input/Examples/CorrectCovmat.txt");
  Covfile.open("../input/Examples/ApproximateCovmat.txt");
  Simfile.open("../input/Examples/Simulations.txt");
  
  
   // For i < 50, the mean is mu1 = 1, for the other it is mu2 = 0
  
  //the data
    for(int i = 0; i < ndat; i++)
    {
      
      if(i < 50)
      {
       Datfile << gsl_ran_gaussian (r,1.0) + 1.0 << " ";
      }
      
      else
      {
        Datfile << gsl_ran_gaussian (r,1.0) + 0.0 << " ";
      }
      
    }
    
    Datfile << endl;
    
    //covmat
    double eps = 1e-2;
    for(int i = 0; i < ndat; i++)
    {
     for(int j = i; j < ndat; j++)
     {
       if(i == j)
       {
         Covfile << i+1 << " " << j+1 << " " << 1.0 + eps*abs(gsl_ran_gaussian (r,1.0)) << endl; 
       }
      
       else
       {
         double c = 0.0 + eps*abs(gsl_ran_gaussian (r,1.0)) ;
         Covfile << i+1 << " " << j+1 << " " << c << endl; 
         Covfile << j+1 << " " << i+1 << " " << c << endl; 
       }
     }
    }
    
  
  //the simulations of the data
  for(int z = 0; z < nsims; z++)
  {
    for(int i = 0; i < ndat; i++)
    {
      
      if(i < 50)
      {
       Simfile << gsl_ran_gaussian (r,1.0) + 1.0 << " ";
      }
      
      else
      {
        Simfile << gsl_ran_gaussian (r,1.0) + 0.0 << " ";
      }
      
    }
    
    Simfile << endl;
  }
  
  Datfile.close();
  Covfile.close();
  Simfile.close(); 
  
  
  /*
   * Step 2: we measure the coverage
   */
  
  
  vector<double> TrueParams(2);
  TrueParams[0] = 1;
  TrueParams[1] = 0; 
  
  
  Example1 E1("../input/Examples/DataVector.txt", "../input/Examples/Simulations.txt", 100, nsims, TrueParams, "Example1");
  //E1.LoadCovmat_KiDSFormat("../input/Examples/CorrectCovmat.txt");
  E1.ComputeCoverage();
  
  
  
}
