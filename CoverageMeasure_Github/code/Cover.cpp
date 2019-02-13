/*
 * Source code to measure coverages of posteriors.
 * If you have faithful simulations and a good approximation
 * for your likelihood function, then this code allows you to
 * check whether your data analysis causes biases. The output
 * files allow you to correct for those biases.
 * 
 * Elena Sellentin, 2018.
 * Sterrewacht
 * Universiteit Leiden
 * The Netherlands
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include "./utils/Matrices.h"
#include "Cover.h"
using namespace std;


/**
 * True params must be the values of the parameters with which the
 * simulations were run.
 * 
 * data_vector_file is supposed to be a single column of numbers,
 * each number one data point.
 * 
 * simulations_file is expected to be a blank-separated file,
 * organized as a matrix: each row one simulation.
 * Nsim rows in total, since there are Nsim simulations.
 * Each line is expected to contain blank-separated the
 * Ndat simulated data points. Meaning each row is Ndat long.
 */
Cover::Cover(string data_vector_file, string simulations_file, int Ndata_intended, int Nsims_intended,vector<double>TrueParams, string suffix_in):Ndat(Ndata_intended), Nsim(Nsims_intended),TrueParamValues(TrueParams), suffix(suffix_in)
{
    /*
     * This constructor does a lot of sizing,
     * since there will be OpenMP routines
     * operating on the vectors later on.
     * This is only possible if there is
     * no resizing, push_back, insert, etc 
     * happening. (Otherwise, the index-lookup
     * is not guaranteed to be flushed in time
     * for OpenMP to do what we want.)
     */
    
    cout << endl << endl << "~~~~~~~~~~~~~~~~~ Coverage Measurement ~~~~~~~~~~~~~~~~" << endl;
    
    OutputFolder(suffix);
    
    data.resize(Ndat);
    simulations.resize(Nsim);
    for(int i = 0; i < Nsim; i++)
    {
        simulations[i].resize(Ndat);
    }
    
    
    ifstream datfile;
    datfile.open(data_vector_file.c_str());
    
    ifstream simfile;
    simfile.open(simulations_file.c_str());
    
    if(!datfile.is_open()) {cout << "Error: Can't find input file for data." << endl; abort();}
    if(!simfile.is_open()) {cout << "Error: Can't find input file for simulations." << endl; abort();}
    
    for(int i = 0; i < Ndat; i++)
    {
        double temp;
        datfile >> temp;
        data[i] = temp;
    }
    
    for(int N = 0; N < Nsim; N++)
    {
        for(int i = 0; i < Ndat; i++)
        {
            double temp;
            simfile >> temp;
            simulations[N][i] = temp;
        }
    }
    
    
    int nop = TrueParamValues.size();
    pof = nop;
    
    //used to track how the best fit scatters around when data are exchanged
    //interesting, but not of central importance
    BestFitParams.resize(Nsim+1);
    for(int i=0; i < Nsim+1; i++)
    {
        BestFitParams[i].resize(nop);
    }
    
    RequestTheseContourLevels();
    int c = contourlevels.size();
    
    
    posteriorheight.resize(Nsim+1);
    for(int i = 0; i < Nsim+1; i++)
    {
      posteriorheight[i].resize(nop);
      for(int p = 0; p < nop; p++)
      {
        posteriorheight[i][p].resize(nop);
        for(int b = p; b < nop; b++)
        {
          //OpenMP will write to these vectors in parallel, so be careful.
            posteriorheight[i][p][b].resize(c);
            for(int cip = 0; cip < c; cip++)
            {
              posteriorheight[i][p][b][c] = 0;
            }
        }
      }
    }
    
    
    covervalues.resize(nop);
    for(int i = 0; i < nop; i++)
    {
        covervalues[i].resize(nop);
        for(int j = i; j < nop; j++)
        {
            //OpenMP will write to these vectors in parallel, so be careful.
            covervalues[i][j].resize(c);
            for(int cip = 0; cip < c; cip++)
            {
              covervalues[i][j][c] = 0;
            }
        }
    }
    
    PosteriorAtTruth.resize(Nsim+1);
    
    C = gsl_matrix_calloc(Ndat,Ndat);
    invC = gsl_matrix_calloc(Ndat,Ndat);
    Corr = gsl_matrix_calloc(Ndat,Ndat);
    
    gsl_matrix_set_identity(C);
    gsl_matrix_set_identity(invC);
    gsl_matrix_set_identity(Corr);
    
    grid_setup = false;
    ran_checks = false;
    make_report = true;
    
}



/**
 * This matrix loader is very specific for the KiDS team format.
 */
void Cover::LoadCovmat_KiDSFormat(string matrixfile)
{

    
    ifstream in;
    in.open(matrixfile.c_str());
    if(!in.is_open())
    {
        cout << "Can't find input file for covariance matrix in KiDS format" << endl;
        abort();
        
    }
    string trash; //KiDS covmats start with a header, which we trash here
    getline(in,trash);
    while(in.good())
    {
        int i,j;
        double temp;
        in >> i; in >> j; in >> temp;
        gsl_matrix_set(C, i-1, j-1, temp); //the minus 1 arises since C counts from zero, but KiDS from 1.
        gsl_matrix_set(C, j-1, i-1, temp);
        
        //cout << i-1 << " " << j-1 << " " << temp << endl;
    }

    for(int i = 0; i < Ndat; i++)
    {
        for(int j = 0; j < Ndat; j++)
        {
            double temp = gsl_matrix_get(C,i,j);
            temp = temp/sqrt( gsl_matrix_get(C,i,i)*gsl_matrix_get(C,j,j));
            gsl_matrix_set(Corr, i,j, temp);
        }
    }
    
    invmatrix(C, Ndat, invC);
    
}


// 

double Cover::Likelihood(vector<double> const something_eg_mean, vector<double> const data_or_sims)
{
    if(something_eg_mean.size() != data.size())
    {
        cout << "This example assumes you use 'something' to be the mean of a Gaussian, and 'data' to be some data vector. So they must have the same dimension. Please derive this class and change the setup, if this is not what you wish to have." << endl;
        abort();
    }
    
    return exp(-0.5* chisqCvecsInvmat(data_or_sims, something_eg_mean, invC) );
    
}



double Cover::Likelihood_on_Grid(int const p1, int const p2, int const gp_1, int const gp_2, int const dat_o_sim)
{
    if( (dat_o_sim < 0 ) || (dat_o_sim > Nsim) )
    {
        cout << "You called the likelihood for something that is neither the true data, nor a simulated data set. Please fix the integers which you feed into Likelihood_on_Grid(int, int, problematic int)." << endl;
        abort();
    }
    
    
    double res = 0;
    if(dat_o_sim == 0)//then take the real data
    {
        res= Likelihood(TheoryAtGridPoints[p1][p2][gp_1][gp_2], data);
    }
    
    if(dat_o_sim > 0)//then take simulated data
    {
        
        res= Likelihood(TheoryAtGridPoints[p1][p2][gp_1][gp_2], simulations[dat_o_sim-1]);
    }
    
    return res;
    
}
  
  
  
double Cover::Prior(vector<double> const parameters)
{
    //flat prior per default
    return 1.0;
    
    //could also be a Gaussian or anything else
    /*
     vector<double> means(TrueParamValues.size());
     means[0] =
     means[1] =
     vector<double> variances(TrueParamValues.size());
     variances[0] =
     variances[1] =
     
     double prior = 1.0;
     for(int i = 0; i < means.size(); i++ )
     {
       prior = prior*(  exp(-0.5*pow(parameters[i] - means[i],2)/variances[i]) / sqrt(variances) ); //ignores 2pi in prior norm
     }
     
     return prior;
     */
     
}

  
  
double Cover::Prior_on_Grid(int const p1, int const p2, int const gp_1, int const gp_2)
{
    vector<double> parameter(2);
    parameter[0] = CoordToPhysical(p1, gp_1);
    parameter[1] = CoordToPhysical(p2, gp_2);
    
    
    return Prior(parameter);
}
  

/*
 * Users, please do not modify this routine. 
 * The functions intended to be modified are marked by "virtual". 
 */
void Cover::Posterior_on_Grid(int const p1, int const p2, int const dat_o_sim)
{
    
    if(p2 < p1)
    {
        cout << "We only provide storage space for posterior values where p1 geq p2, due to symmetry reasons. Please interchange your arguments." << endl;
        abort();
    }
    
    if( (dat_o_sim < 0 ) || (dat_o_sim > Nsim) )
    {
        cout << "You called the likelihood for something that is neither the true data, nor a simulated data set. Please fix the integers which you feed into Likelihood_on_Grid(int, int, problematic int)." << endl;
        abort();
    }
    
    for(int gp_1=0; gp_1 < inc[p1]; gp_1++)
    {
        for(int gp_2=0; gp_2 < inc[p2]; gp_2++)
        {
            
            AllPosteriors[dat_o_sim][p1][p2][gp_1][gp_2] = Likelihood_on_Grid(p1, p2, gp_1,gp_2,dat_o_sim)*Prior_on_Grid(p1,p2,gp_1,gp_2);
            
            if( iskrampf(AllPosteriors[dat_o_sim][p1][p2][gp_1][gp_2], __func__) )
            {
              cout << "Nonsensical posterior for dat_o_sim " << dat_o_sim << " for parameter pair " << p1 << "," << p2 << " for grid point " << gp_1 << "," << gp_2 << endl;
              abort();
            }
            
            if( AllPosteriors[dat_o_sim][p1][p2][gp_1][gp_2] < 0 )
            {
              cout << "Negative posterior for dat_o_sim " << dat_o_sim << " for parameter pair " << p1 << "," << p2 << " for grid point " << gp_1 << "," << gp_2 << endl;
            }
    
        }
    }
}



double Cover::Posterior_at_TrueParams(int const dat_o_sim)
{
    
   double res = 0;
   if(dat_o_sim == 0)
   { 
       res = Likelihood(TheoryAtTrueParamPoint, data)*Prior(TrueParamValues);
   }
   
   if(dat_o_sim > 0)
   { 
       res = Likelihood(TheoryAtTrueParamPoint, simulations[dat_o_sim-1])*Prior(TrueParamValues);
   }
   return res;
}
  
  
  
void Cover::AssignTheoryToGridFromFiles(string TreePrefix, vector<vector<string> > Foldertree, vector<string> fileprefixes, string TrueParamFile)
{
   //Do not edit the source code here; derive and edit there.
}  
  
  
void Cover::AssignTheoryToGridFromFunction(int p1, int p2, int gp_1, int gp_2)
{
   //Do not edit the source code here; derive and edit there.
}
  
void Cover::SetupGrid(vector<double> lowerBin, vector<double> upperBin, vector<int> incin)
{
    
    
    if(lowerBin.size() != upperBin.size())
    {
        cout << "Size mismatch in SetupGrid: your vectors have different dimensions." << endl;
        abort();
    }
    
    if(lowerBin.size() != incin.size())
    {
        cout << "Size mismatch in SetupGrid: your vectors have different dimensions." << endl;
        abort();
    }
    
    lowerB.resize(0);
    upperB.resize(0);
    inc.resize(0);
    
    for(int i = 0; i < lowerBin.size(); i++)
    {
        lowerB.push_back(lowerBin[i]);
        upperB.push_back(upperBin[i]);
        inc.push_back(incin[i]);
    }
    
    for(int i = 0; i < incin.size(); i++)
    {
      delta.push_back( (upperB[i] - lowerB[i])/( (double)inc[i]) );
    }
    
    
    
    
    
    /*
     * Now we provide storage space for all the 
     * posteriors which we need to compute. This
     * can be huge! If we didn't do this, we'd
     * need to always recompute things and would
     * instead be highly CPU intensive.
     * 
     * Note that OpenMP writes to different 
     * elements of this array at the same time.
     */
    
    TheoryAtGridPoints.resize(pof);
    for(int p1 = 0; p1 < pof; p1++)
    {
      TheoryAtGridPoints[p1].resize(pof);
      
      for(int p2 = p1; p2 < pof; p2++)
      {
        TheoryAtGridPoints[p1][p2].resize(inc[p1]);
        for(int i = 0; i < inc[p1]; i++)
        {
          TheoryAtGridPoints[p1][p2][i].resize(inc[p2]);
        }
      }
      
    }
    
    AllPosteriors.resize(Nsim+1); //because index i = 0 will store the posterior of the real data
    for(int i=0; i <= Nsim; i++)
    {
        AllPosteriors[i].resize(pof);
        for(int p1=0;p1<pof;p1++)
        {
            AllPosteriors[i][p1].resize(pof);
            for(int p2 = p1; p2 < pof; p2++) //notice that we onlz provide space for p2 \geq p1 !
            {
                AllPosteriors[i][p1][p2].resize(inc[p1]);
                for(int g=0; g < inc[p1]; g++)
                {
                    AllPosteriors[i][p1][p2][g].resize(inc[p2]);
                }
            }
        }
        
    }
    
    grid_setup = true;
}
  
  
  
void Cover::ComputeCoverage()
{
  
    
    if(grid_setup == false)
    {
        cout << "Please provide upper and lower boundaries for the posterior grid, by calling the functing SetupGrid(vec, vec, vec)." << endl;
        abort();
    }
    if(ran_checks == false){ Checks();}
    if(ran_checks == false)
    {
        cout << "The function Checks() usually sets the bool ran_checks to true, after having successfully completed. It does not do so anymore, so unintended changes to the code must have occured. Please fix this." << endl;
        abort();
    }
    
   
    Checks();
    
    /*
     * First compute the posterior values for the 
     * true parameters. If these Ps are larger than
     * the height where the confidence contours cut,
     * then the coverage counter is increased.
     * 
     * FAILURE MODE: if users normalize the posterior
     * on the grids, but forget to then divide the 
     * posteriorheight at the true params by the same
     * number. We do not normalize by default.
     * We hope that by having made the function Posterior_on_Grid
     * non-overwritable will secure anyone against this trap.
     */
    

    
    for(int N = 0; N < Nsim+1; N++)
    {
        PosteriorAtTruth[N] = Posterior_at_TrueParams(N);
        if(PosteriorAtTruth[N] < 0)
        {
          cout << "Your code produces negative posterior values." << endl;
          abort();
        }
    }
    
    //the number of params to be varied
    int nop = TrueParamValues.size();
    
    
    cout << "Starting to measure coverage for parameters..." << endl;
    
    //these two loops run over the corner plot
    //they select pairs of parameters
    for(int p1=0; p1<nop; p1++)
    {
        for(int p2=p1; p2<nop; p2++)
        {
          cout << "..." << p1 << " " << p2 << endl;
            
            /**
             * ATTENTION: might be that I'll need to provide a function-internal version of 
             *  my current member variable AllPosteriors here. OpenMP standard does
             *  not by default allow parallel computation on a member variable (because
             *  that member variable might not have been instantiated yet). OpenMP outlines
             *  parallel regions...
             */
            
            if(AllPosteriors[0][p1][p2].size() == 0)
            {
                cout << "Something went wrong in how the vector AllPosteriors is sized. Please correct the problem." << endl;
                abort();
            }
            
            if(p1==p2) 
            {
                continue; // because I have not yet coded up the 1d case
            }
            
            
            if(p1!=p2)
            {
                
               /*
                * Now we have selected a pair of two distinct parameters.
                */
                TwoDimensionalCase(p1,p2);
                string covername = suffix+"/CoverageTable"+suffix;
                CoverageTableToFile(p1, p2, covername);
            }//p1 \neq p2
            
            
        }//p2-loop
    }//p1-loop
    
    make_report = true;
}



void Cover::GridToVector(int const p1, int const p2, int dat_o_sim, vector<double> & output)
{
    
    output.resize(0);
    int xdim = AllPosteriors[dat_o_sim][p1][p2].size();
    int ydim = AllPosteriors[dat_o_sim][p1][p2][0].size();
    
    for(int i = 0; i < xdim; i++)
    {
        for(int j = 0; j < ydim; j++)
        {
            output.push_back(AllPosteriors[dat_o_sim][p1][p2][i][j]  );
        }
    }
}



void Cover::TwoDimensionalCase(int const p1, int const p2)
{
    
    /*
     * For fixed parameter pair p1, p2, loop over all simulations,
     * and for N = 0, the real data. The loop computes posteriors.
     * 
     * When all posteriors were computed,
     * compute all their confidence contours.
     */
    for(int N = 0; N < Nsim+1; N++)
    {
        Posterior_on_Grid(p1,p2, N);
        
        vector<double>testprob;
        
        GridToVector(p1, p2, N, testprob);
        
        OutputFolder(suffix+"/Posteriors");
        
        string name = suffix+"/Posteriors/Contours" + strint(p1)+"_"+strint(p2)+"_"+strint(N);

        if(N == 0)//then it is the posterior of the real data, which we definitely need
        { 
          ConfidenceContours(testprob, posteriorheight[N][p1][p2],name);
          PosteriorsToFile(p1, p2, N, "Post");
        }
        
        
        else if(N % 50 == 0)//a few other posteriors, to look around
        { 
          ConfidenceContours(testprob, posteriorheight[N][p1][p2],name);
          PosteriorsToFile(p1, p2, N, "Post");
        }
        
        else
        {
          ConfidenceContours(testprob, posteriorheight[N][p1][p2],"");
        }
    }
    
    
    /*
     * Now we have where the default analysis puts the confidence contours
     * for all posteriors. Now we measure the coverage. We measure it
     * from simulations ONLY. (Later, we APPLY it to the data posterior).
     * Hence the loop range now starts at 1, because zero is for the real data.
     */
    
    for(int N = 1; N < Nsim+1; N++) //starting at 1 to exclude the data
    {
        for(int c = 0; c < contourlevels.size(); c++)
        {
            //c is the index determining the contour level; Careful with OpenMP!
            if(PosteriorAtTruth[N] >= posteriorheight[N][p1][p2][c] ) //then contour contains true point
            { 
                covervalues[p1][p2][c] = covervalues[p1][p2][c] + 1; 
            }
        }
    }
    
    //divide by Nsim to get percentages
    for(int c = 0; c < contourlevels.size(); c++)
    {
        covervalues[p1][p2][c] = (covervalues[p1][p2][c]/(double)(Nsim) );
    }
    
}
    
void Cover::OneDimensionalCase(int const p1)
{
   //TODO 
}
    
    
void Cover::RequestTheseContourLevels()
{
    /*
     * We specify default contour levels,
     * of which the coverage is computed.
     * The member variable 'contourlevels'
     * stores these.
     * 
     * This is not a virtual function,
     * since we only want users who wish to extend
     * the method to modify the requested contour
     * levels. Exchanging physical applications
     * is independent of this function.
     */
    
    contourlevels.resize(0);
    
    for(int i=0; i<100; i++)
    {
        contourlevels.push_back(i);
    }
    
    
    /*
     * Now investigate the outskirts of
     * the posterior a bit more precisely.
     * Note that this requires a well resolved
     * parameter grid. Otherwise, the outskirts
     * are badly resolved.
     */
    for(int i=0; i<10; i++)
    {
        contourlevels.push_back(90.5 + i);
    }
    
    for(int i=0; i<5; i++)
    {
        contourlevels.push_back(95.25 + i);
        contourlevels.push_back(95.75 + i);
    }
    
    sort(contourlevels.begin(),contourlevels.end());
    
    for(int i = 0; i < contourlevels.size(); i++)
    {
        if(contourlevels[i] > 99.99)
        {
            cout << "Warning: you should not request such high contour levels as you did. It is highly unlikely that the simulations and the grid resolution are able to reliably resolve the 100% confidence contour. Should you have a specific setup where this is possible (like compact support for all distributions), then please delete this check." << endl;
            abort();
        }
        
        if(i < contourlevels.size()-2) //otherwise i+1 goes beyond the end of the vector
        {
            if(contourlevels[i] > contourlevels[i+1])
            {
                cout << "Your requested contour levels are not in monotonically increasing order. Please fix this, the code depends on it." << endl;
                abort();
            }
        }
    }
    
    
    /*
     * Given Nsim and the above determined contour levels,
     * we can already compute the standard deviation of the 
     * coverage: the standard deviation is the sqrt of
     * the variance of a Binomial process, see paper,
     * where "fraction" is "alpha" in the paper.
     */
    for(int i = 0; i < contourlevels.size(); i++)
    {
        double fraction = contourlevels[i]/100.0; // fraction in [0,1]
        
        covervalues_stdv.push_back( sqrt( fraction*(1.-fraction)/(double)(Nsim)    )  );
    }
}
    

void Cover::ConfidenceContours(vector<double> const testprob, vector<double> & heights, string saveto)
{
    /*
     * Copy the incoming probability vector, 
     * since the copy will be reshuffled.
     */
    vector<double> pdf(testprob);
    
    
    heights.resize(contourlevels.size());

    double sum = 0;
    double sumsig = 0;
    double contsig = 0;


    //normalizing and shuffling
    for(int i= 0; i < pdf.size(); i++)
    { 
        if(pdf[i] == 0)
        {
          cout << "Warning: posterior likelihood smaller than eps of C++ (meaning it is treated as zero, and this will affect the confidence levels)." << endl;
        }
        sum += pdf[i]; 
    }
    sort(pdf.begin(), pdf.end());
    
    
    int start = 0;

    if(contourlevels[start] == 0)
    {
        heights[0] = pdf[pdf.size() -1]; //The zero-percent confidence contour is where the peak is
        start++;
    }
    
    //Now start at the zeroth or first confidence contour
    //depending on what start is by now
    double conflev = contourlevels[start];
    
    /*.size-1 because arrays count from zero*/
    for(int i = pdf.size() - 1; i >= 0; i--)
    {
      //cout << i << endl;
      
      
        sumsig += pdf[i];
        contsig = pdf[i];

        if (sumsig > sum * conflev/100.0)
        { 
            heights[start] = contsig;
            start++;
            
            if(start == contourlevels.size())
            {
                //if, then there are no more levels to compute, hence break.
                break;
            }
            
            
            if(start != contourlevels.size())
            {
              
                //then there are further levels so
                //go to the next contour
                conflev=contourlevels[start];
            }
        }
    }

    ofstream of;
    of.open(saveto.c_str());
    of << "# contour level      posterior height" << endl;
    for(int l = 0; l < contourlevels.size(); l++)
    {
        of << contourlevels[l] << "        " << heights[l] << endl;
    }
    
    of.close();

}





    
void Cover::Checks()
{
 
    int pdim = TrueParamValues.size();
    
    if(lowerB.size() != pdim)
    {
        cout << "Error: There is not a lower boundary for all parameter values. Please change TrueParams (constructor) or lowerB (SetupGrid)." << endl;
        
        abort();
    }
    
    if(upperB.size() != pdim)
    {
        cout << "Error: There is not an upper boundary for all parameter values. Please change TrueParams (constructor) or upperB (SetupGrid)." << endl;
        
        abort();
    }
    
    
    for(int i = 0; i < pdim; i++)
    {
        if(lowerB[i] > TrueParamValues[i])
        {
            cout << "Error: for parameter " << i << " the lower boundary excludes the true parameter value." << endl;
            abort();
        }
        
        if(upperB[i] < TrueParamValues[i])
        {
            cout << "Error: for parameter " << i << " the upper boundary excludes the true parameter value." << endl;
            abort();
        }
    }
    
    
    int nof_zeroes = 0;
    
    for(int N = 0; N < Nsim; N++)
    {
        //Using here that c++ has preimplemented routines to compare
        //doubles to zero, even though it looks like a bug par excellence
        
        for(int i = 0; i < Ndat; i++)
        {
          if(simulations[N][i] == 0){nof_zeroes++;}
        }
    }
    
    if(nof_zeroes > 0.02*Nsim*Ndat)
    {
        cout << "There are surprisingly many zeroes in your simulations. Are you sure you read in the simulations correctly? If sure, then take out this check and recompile. This is File Cover.cpp, function Checks()." << endl;
        abort();
    }
    
    int nof_dat_zeroes = 0;
    for(int i = 0; i < Ndat; i++)
    {
        if(data[i] == 0){ nof_dat_zeroes++; }
    }
    
    if(nof_dat_zeroes > 0.02*Ndat)
    {
        cout << "There are surprisingly many zeroes in your data. Are you sure you read in the data correctly? If sure, then take out this check and recompile. This is File Cover.cpp, function Checks()." << endl;
        abort();
    }
    
    
    if(contourlevels.size() == 0)
    {
        cout << "Cannot compute coverage: no contourlevels were given. Please run the function 'RequestTheseContourLevels()' first." << endl;
        abort();
    }
    
    
    for(int i = 0; i < contourlevels.size(); i++)
    {
        if(contourlevels[i] > 99.99)
        {
            cout << "Warning: you should not request such high contour levels as you did. It is highly unlikely that the simulations and the grid resolution are able to reliably resolve the 100% confidence contour. Should you have a specific setup where this is possible (like compact support for all distributions), then please delete this check." << endl;
            abort();
        }
        
        if(i < contourlevels.size()-2) //otherwise i+1 goes beyond the end of the vector
        {
            if(contourlevels[i] > contourlevels[i+1])
            {
                cout << "Your requested contour levels are not in monotonically increasing order. Please fix this, the code depends on it." << endl;
                abort();
            }
        }
    }
    
    
    if(TheoryAtTrueParamPoint.size() == 0)
    {
        cout << "You don't seem to have provided any theory input for the true parameters of the simulations. The vector AtTrueParamPoint has size zero." << endl;
        abort();
    }
    
    
    
    ran_checks = true;
    
}


void Cover::FreeMemory()
{
    if(invC != NULL)
    {
      gsl_matrix_free(invC);
    }
    
    if(C != NULL)
    {
      gsl_matrix_free(C);
    }
    
    if(Corr != NULL)
    {
      gsl_matrix_free(Corr);
    }
    
}



void Cover::PosteriorsToFile(int p1, int p2, int dat_o_sim, string  name)
{
  
  OutputFolder(suffix+"/Posteriors");

  string filename = suffix+"/Posteriors/"+ name + strint(p1) + "_" + strint(p2) + "_" + strint(dat_o_sim);
  
  ofstream posts;
  posts.open(filename.c_str());
  for(int i = 0; i < inc[p1]; i++)
  {
    for(int j = 0; j < inc[p2]; j++)
    {
      
      posts << CoordToPhysical(p1,i) << " " << CoordToPhysical(p2,j) << " " << AllPosteriors[dat_o_sim][p1][p2][i][j] << endl;
    }
    
    posts << endl;
  }
  
  posts.close();
  
}




/*Have to guard this against going out of grid boundaries*/
double Cover::CoordToPhysical(int const p, int const gp)
{
    if(gp > inc[p])
    {
        cout << "Parameter point beyond precomputed grid requested, for parameter " << p << " at unavailable increment " << gp << endl;
        abort();
    }
    else
    return lowerB[p] + gp*delta[p];
}


void Cover::CoverageTableToFile(int const p1, int const p2, string filename)
{
    bool mismatch = false;
    
    if(contourlevels.size() != posteriorheight[1][p1][p2].size()){mismatch = true;}
    if(contourlevels.size() != covervalues[p1][p2].size()){mismatch = true;}
    if(contourlevels.size() != covervalues_stdv.size()){mismatch = true;}
    
    if(mismatch)
    {
        cout << "Error: The coverages and the recalibration could not be computed, since the lookup table has nonesensical dimensions. Please correct the lookup table: all vectors must have the same dimension, and for each index i, the same level needs to be targeted." << endl;
        
        abort();
    }
    
    
    ofstream covout;
    covout.open(filename.c_str());
    
    covout << "# Coverage recalibration table. The last column are conservative estimates of the uncertainties." << endl;
    
    covout << "# The first column denotes the level of the confidence or credibility contour \n as it currently arises from the statistical method (e.g. your likelihood approximation) which you have put in. " << endl;
    
    covout << "# The second column stores the coverage probability of your confidence contours. If the numbers in the first and second line do not agree, then your current parameter plots cannot naively be interpreted as 'the 68% confidence contour contains the true value of the parameters 68% of the time, if my model is correct.' Instead, the contour will contain the true parameters more often or less often that stated by the label of your contour." << endl;
    
    covout << "# Confidence or credibility contours cut the posterior at constant height. The third column stores the posterior height at the confidence levels given in the first column. To be precise, the third column stores P(theta|data), being the posterior for the REAL data, not the simulated ones. If you wish to plot confidence contours with correct coverage, you have to search for the coverage value of your choice in column 2, and cut the posterior at the value stated in column 3. This is what the automatically produced output plots do." << endl;
    
    
    covout << endl << endl;
    covout << "#contour-level [0,1]    coverage probability [0,1]         posterior height            std_dev of column 2" <<endl;
    
    
    /*
     * TODO: posteriorheight[0] should be the data posterior. Is it?
     */
    for(int i=0; i < contourlevels.size(); i++)
    {
       covout << setprecision(5) << contourlevels[i]/100. <<  "              " << covervalues[p1][p2][i] << "       ";
       covout << setprecision(5) << posteriorheight[0][p1][p2][i];
       covout << setprecision(3) << "       " << covervalues_stdv[i] << endl;
    }
    
    covout.close();
    
    
    
    ofstream Gnup;
    string Gnupname = suffix+"/CoverTablePlotter_"+suffix+".gp";
    Gnup.open(Gnupname.c_str());
    
    Gnup << "set term postscript color solid enhanced eps" << endl;
    Gnup << "set out \""+suffix+"/CoverTable"<< suffix<<".eps\"" << endl;
    Gnup << "set grid; set size 0.65,0.7; set ylabel 'Coverage probability'; set xlabel 'Confidence or credibility contour label'" << endl;
    
    
    Gnup << "plot '" << filename << "'  u 1:2:4 w err pt 7 ps 0.7 lc rgb 'blue' title 'measured coverage', x lc rgb 'navy' title 'Perfect coverage'" << endl;
    Gnup << "set out " << endl;
    
    Gnup.close();
    
    string command = "gnuplot " + Gnupname;
    system(command.c_str());
    
    
    ApplyCoverageToDataPosterior();
}




void Cover::ApplyCoverageToDataPosterior()
{
    
  for(int i = 0; i < pof; i++)
  {
    for(int j = i+1; j < pof; j++)
    {
      
      int xcount = AllPosteriors[0][i][j].size();
      int ycount = AllPosteriors[0][i][j][0].size();
      string datafile = suffix+"/Posteriors/Post"+ strint(i) + "_" + strint(j) + "_" + strint(0);
      string exname = suffix+"/Python"+strint(i)+"_"+strint(j);
      
      pythoncontours( xcount, ycount, exname, datafile, i, j);
      
       
    }
  }
    
}






void Cover::pythoncontours( int xcount, int ycount, string executablename, string datafile, int p1, int p2)
{

    double contsig1 = 0;
    double contsig2 = 0;
    double contsig3 = 0;
    
    double uncorr1 = 0;
    double uncorr2 = 0;
    double uncorr3 = 0;


    for(int i=0; i < contourlevels.size(); i++)
    {
        contsig1 = posteriorheight[0][p1][p2][i];

        if (covervalues[p1][p2][i] > 0.68) break;
    }
    
    
    for(int i=0; i < contourlevels.size(); i++)
    {
        contsig2 = posteriorheight[0][p1][p2][i];

        if (covervalues[p1][p2][i] > 0.90) break; 
    }
    
    
    for(int i=0; i < contourlevels.size(); i++)
    {
        
        contsig3 = posteriorheight[0][p1][p2][i];

        if (covervalues[p1][p2][i] > 0.95) break; 
    }



    for(int i=0; i < contourlevels.size(); i++)
    {
        uncorr1 = posteriorheight[0][p1][p2][i];

        if (contourlevels[i]> 68) break; 
    }
    
    for(int i=0; i < contourlevels.size(); i++)
    {
        uncorr2 = posteriorheight[0][p1][p2][i];

        if (contourlevels[i]> 90) break; 
    }
    
    
    for(int i=0; i < contourlevels.size(); i++)
    {
        uncorr3 = posteriorheight[0][p1][p2][i];

        if (contourlevels[i]> 95) break; 
    }
    
  

    string xaxis = strint(p1);
    string yaxis = strint(p2);

    FILE * pythonex;
    pythonex = fopen(executablename.c_str(), "w");

    fprintf(pythonex, "import numpy as np\nimport matplotlib.pyplot as plt\nfrom matplotlib import rc\n\n" );
    
    fprintf(pythonex, "from matplotlib import rcParams\nrcParams.update({'figure.autolayout': True})\n");
    
    fprintf(pythonex, "import matplotlib.patches as mpatches\n");
    fprintf(pythonex, "font ={'size': 15}\nplt.rc('font',**font)\nplt.rc('text',usetex=True)\norigin = 'lower'\n\n");

    fprintf(pythonex, "N2x=%.d\nN2y=%.d\n",xcount,ycount);
    fprintf(pythonex, "x2,y2,z2 = np.genfromtxt(r'%s', unpack=True)\n",datafile.c_str());

    fprintf(pythonex, "xi2 = np.reshape(x2,(-1,N2x))\nyi2 = np.reshape(y2,(-1,N2y))\nzi2 = np.reshape(z2,(-1,N2x))\n");

    //Some matplotlib installs want the 3 contour levels to be given in reverse order
    fprintf(pythonex, "levels2 = [%.3e,%.3e,%.3e]\n",contsig3, contsig2, contsig1);
    fprintf(pythonex, "CS2 = plt.contour(xi2, yi2, zi2, levels2, colors = ('mediumblue','blue'),linewidths = (2),origin = origin)\n");
    
    fprintf(pythonex, "red_patch = mpatches.Patch(color='mediumblue', label='coverage corrected')\n" );
    fprintf(pythonex, "plt.legend(handles=[red_patch])\n\n");
    
    
    
    //Some matplotlib installs want the 3 contour levels to be given in reverse order
    fprintf(pythonex, "levelsun = [%.3e,%.3e,%.3e]\n",uncorr3, uncorr2, uncorr1);
    fprintf(pythonex, "CS2 = plt.contour(xi2, yi2, zi2, levelsun, colors = ('gold','goldenrod'),linewidths = (1.8), alpha = 0.9, origin = origin)\n");
    
    
    fprintf(pythonex, "plt.xlabel(r'%s')\nplt.ylabel(r'%s')\n",xaxis.c_str(),yaxis.c_str());
    fprintf(pythonex, "\n");

    fprintf(pythonex, "plt.plot([%.3e], [%.3e], 'bo')\n",TrueParamValues[p1],TrueParamValues[p2]);
    
   
    
    string PDFName =  suffix+"/"+strint(p1)+"_"+strint(p2)+"_covercorrected" +(string)".pdf";
    fprintf(pythonex, "plt.savefig('%s')\n",PDFName.c_str());

    fclose(pythonex);
    
    cout << endl << endl << "Any complaints now following are due to python; please check your local installation, if any arise." <<endl;
    
    string command = "python " + executablename;
    system(command.c_str() );
    cout << endl<< "Coverage corrected posterior written to " << datafile << ".png" << endl;

}







void Cover::EnglishSummary(string filename)
{
    string Summary = "~~~~~~~~ SUMMARY ~~~~~~~~~\n";
    string NotConservative = "Your statistical analysis is *not* conservative.\nThe confidence contours undercover, meaning the true point lies mostly outside of them. The confidence contours are therefore very likely too small. If you used a purposefully informative prior, then this is fine. In all other cases, it is an indication that your statistical analysis must somewhere contain an unrecognized source of systematics.";
    
    string Conservative = "Your statistical analysis is likely on the conservative side.\nThe confidence contours overcover, meaning they are so wide that they contain the true point surprisingly often. If your simulations are fully representative of your real data, then you likely overestimate your data variability during your analysis.";
    
    
    
}



Cover::~Cover()
{
    
    if(make_report){ EnglishSummary(); }
    
    FreeMemory();
    
    cout << endl << endl;
}

