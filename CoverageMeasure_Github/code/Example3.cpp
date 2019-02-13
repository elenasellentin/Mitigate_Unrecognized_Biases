#include "Cover.h"

#include<vector>
#include<cmath>
#include<iostream>
using namespace std;

/*
 * We derive from the base class, and implement
 * an example where an informative prior was
 * hidden in the analysis, even though an uninformative
 * prior was sought for.
 */




class Example3: public Cover
{
  
private:
protected:
public:
  
    Example3(string data_vector_file, string simulations_file, int Ndata_intended, int Nsims_intended,vector<double>TrueParams, string suffix_in):Cover(data_vector_file,simulations_file,Ndata_intended, Nsims_intended,TrueParams, suffix_in)
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
    
  
  //The unforeseen informative prior
  double Prior(vector<double> const parameters)
  {

      vector<double> means(2);
      means[0] = 1.45;
      means[1] = 1.35;
      vector<double> variances(2);
      variances[0] = 0.1;
      variances[1] = 0.1;
     
      double pri = 1.0;
      for(int i = 0; i < means.size(); i++ )
      {
        pri = pri*(  exp(-0.5*pow(parameters[i] - means[i],2)/variances[i]) / sqrt(variances[i]) ); //ignores 2pi in prior norm
      }
     
      return pri;
      
    
   
   //return 1.0;
  
  }
};




int main()
{

  
  int ndat = 100; //100-dimensional data vector
  int nsims = 2000; //thousand 100 dimensional simulations
  
  vector<double> TrueParams(2);
  TrueParams[0] = 1;
  TrueParams[1] = 0; 
  
  
  Example3 E3("../input/Examples/DataVector.txt", "../input/Examples/Simulations.txt", 100, nsims, TrueParams, "Example3");
  E3.ComputeCoverage();
  
  
  
  
  
  
  
}
