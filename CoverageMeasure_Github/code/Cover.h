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

#ifndef LNA_COVER_H
#define LNA_COVER_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <gsl/gsl_matrix.h>
using namespace std;


class Cover
{
 
private:

     //Set to false in constructor; set to true in SetupGrid().
     bool grid_setup;
     
     bool ran_checks;
     
     bool make_report;
     
     
     /**
      * This stores the posterior values for all cross sections through
      * the posterior that are computed. 
      * The first index selects whether it is the posterior of the real
      * data, or a simulation. AllPosteriors[0][][][][] is for the data,
      * meaning Posterior(p1,p2|data).
      * AllPosteriors[i][][][][] for i in [1,Nsim] is for
      * simulations,meaning Posterior(p1,p2|simulation i).
      * The second and third index select parameters, with convention
      * that p2 \geq p1: Allposteriors[i][p1][p2][][].
      * For p2 < p1, no memory is provided. The last two indices select
      * the grid points for p1 and p2. 
      * 
      * This vector is sized in function SetupGrid()
      */
     vector<vector<vector<vector<vector<double > > > > > AllPosteriors;
     


  
    
protected:
public:
    
      /**
     * Auxilliary matrices, in case a Gaussian likelihood
     * is used.
     * C: the covmat. invC: the inverse covmat.
     * Corr: the correlation matrix.
     * 
     * Important: OpenMP knows that these exist.
     * If you provide other elements to your
     * likelihood, OpenMP might need to be told so.
     * See the functions for likelihoods, priors
     * and posteriors.
     */
     gsl_matrix* C;
     gsl_matrix* invC;
     gsl_matrix* Corr;
  
  
    /* INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    
     string suffix;
     /**
     * At each grid point of the parameter grid,
     * one typically needs to compute sth very costly,
     * which then goes into the likelihood.
     * For a Gaussian likelihood, this costly item
     * is the mean of the data. We hence here provide
     * storage room for those means, such that they do
     * not need to be recomputed all the time.
     * TheoryAtGridPoints[p1][p2][a][b][c] will be the mean of the 
     * cth data point, evaluated for the ath value of
     * the first parameter p1, and the bth value for the 
     * second parameter p2. Hence five indices in total:
     * the first two to select the parameters (Om, w, wa, etc)
     * the second two to select a grid point. The last
     * two select the element of whichever vector is attached to that
     * gridpoint.
     */
     vector<vector<vector<vector<vector<double > > > > >TheoryAtGridPoints;
     
     /*
      * It might be that the true parameters (those of the simulations)
      * are determined by other scientific groups than those who provide
      * the grid. It will then be likely that the true parameters happen
      * to not be a grid point. Therefore we here provide this auxilliary
      * vector, where expensive elements can be stored, which the 
      * likelihood needs. Mind you: if you provide Nsim simulations,
      * then this vector will be used Nsim+1 times, so indeed: store
      * whatever you can.
      */
     vector<double>TheoryAtTrueParamPoint;
    
    
    /*
     * The true parameters: those are
     * the input parameters of the simulations.
     * In the paper indicated by little crosses.
     */
    vector<double> TrueParamValues;
    
    
    /**
     * The contour levels, as you request them.
     * Per default, 120 contour levels are computed:
     * contour levels of zero till 90, and then in steps of 0.5
     * since the outer regions of posterios are often unstable.
     * 
     * It is not recommended to go beyond 150 contour levels, say.
     */
    vector<double> contourlevels; // vector of contour levels, e.g. 95%, 99% etc
    
    /* conservative uncertainty estimate of the 95%'s coverage;
     * only depends on contour level and Nsim
     */
    vector<double> covervalues_stdv; 
    
    /*
     * Upper and lower Boundaries for 
     * the parameter grid.
     * inc is the number of increments in the grid.
     * (Meaning inc is the "number of cells" per axis)
     */
    vector<double> lowerB;
    vector<double> upperB;
    vector<int> inc;
    vector<double> delta;
    
    /**
     * The simulations are stored as a vector of vectors.
     * The first index selects the simulation. If there
     * are Nsim sims, then the first index runs in [0,Nsim-1].
     * The second index runs over the data vector. If there
     * are Ndat datapoints, then the second index runs in
     * [0,Ndat-1]
     */
    vector<vector<double> > simulations;
    
    /**
     * The actual data vector: the one from the sky,
     * or from your lab.
     */
    vector<double> data;
    
    int Nsim, Ndat; //Number of sims, number of data points per sim, which is the same as for the real data
    
    int pof;
    
    
    /* OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    /**
     * output: The coverage probabilities of the
     * (Bayesian or non-Bayesian) credibility regions
     * or confidence regions; naming convention
     * depends on statistical school of thought.
     * _stdv indicates the standard deviation.
     * 
     * The below 4 vectors must has the same
     * dimension, as their elements of same index
     * refer to the elements with the same index
     * in the other vectors.
     * 
     * The have the same length as the vector contourlevels, 
     * see above under section INPUT.
     */
    
    // for dat_o_sim, params p1,p2, the height of the DATA posterior where x% levels cuts
    //posteriorheight[0] is for the real data. posteriorheight[i>0] is for sims.
    vector<vector<vector<vector<double> > > >posteriorheight; 
    
    // for params p1,p2, the coverage (in percentages) of the x% level
    vector<vector<vector<double> > > covervalues;  
    
    /**
     * Not of particular importance and unused by the method,
     * but as a sideeffect this handy diagnostic is produced: 
     * It is a Nsim+1 times Nparams dimensional matrix,
     * storing the bestfitting parameter values per sim and
     * for the real data as well. One can then plot how the best
     * fit scatters around. 
     * BestFitParams[0] stores the best itting params for the 
     * actual data. BestFitparams[1] till [Nsim+1]
     * store it for the simulations.
     */
    vector<vector<double> > BestFitParams;
    
    vector<double> PosteriorAtTruth;
    
    
    /* MEMBER FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    
    /* STEP 1: Provide data and simulations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /**
     * Constructor: provide the filenames of the files which
     * 1) store your data vector
     * 2) store all your simulations
     * 3) provide the number of datapoints
     * 4) provide the number of simulations
     * 5) provide a vector which stores the true parameter values
     */
    Cover(string data_vector_file, string simulations_file, int Ndata_intended, int Nsims_intended, vector<double>TrueParams, string suffix_in);
    
    
    
    

    /* STEP 2: Set up your statistical treatment/ your likelihood */

    /*
     * If you have a Gaussian likelihood (default), provide the covmat,
     * here in KiDS-style ordering.
     */
    void LoadCovmat_KiDSFormat(string matrixfile);
    
    /**
     * gp abbreviates "Grid point", _1 and _2 are the axes.
     * The three functions compute the three densities
     * for all grid points as they were set up by SetupGrid().
     * 
     * The logic of the code is to precompute everything on a grid
     * and store the expensive parts in the auxilliary member variables
     * of this class.
     */
    
    
    /**
     * Virtual functions are intended to be overwritten,
     * if the user wishes to exchange the physical application,
     * but keep the statistics the same.
     * 
     * The vector 'something' can be anything else than the data,
     * which your likelihood needs. Per default, it is programmed
     * to be the mean of a Gaussian.
     */
    //Careful: This function will go through OpenMP, do not introduce data races.
    virtual double Likelihood(vector<double> const something_eg_mean, vector<double> const data_or_sims); 

    
    virtual double Prior(vector<double> const parameters);
    
    /**
     * Derive the class and overwrite this function, if the default 
     * implementation does not match your physical application.
     * Per default, the prior is unnormalized and flat. It simply returns 1.
     */
    //Careful: This function will go through OpenMP, do not introduce data races.
    virtual double Prior_on_Grid(int const p1, int const p2, int const gp_1, int const gp_2);

    /**
     * If dat_o_sim is 0, then the real data are used.
     * If dat_o_sim is [1,Nsim], then simulation (dat_o_sim-1) is used.
     * The -1 arises since data and simulations are stored seperately.
     * In the other cases, an error is thrown.
     * 
     * These functions better be left unchanged by users,
     * except for if you wish to extend the method.
     */
    //Careful: This function will go through OpenMP, do not introduce data races.
    virtual double Likelihood_on_Grid(int const p1, int const p2, int const gp_1, int const gp_2, int const dat_o_sim);
    

    

    
    /* STEP 3: Provide theory input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    
    /**
     * Determines how finely spaced the grid will be. This is non-virtual on purpose!
     * If the grid is not regularly spaced, then volume correction factors would need to
     * go into the confidence contour routines (which they currently don't).
     */
    void SetupGrid(vector<double> lowerBin, vector<double> upperBin, vector<int> spacingin);
    
    /**
     * This function fills the vectors 
     * AtGridPoints and AtTrueParamPoint.
     * Per default, it loads from data files
     * what needs to be stored.
     * 
     * It assumes that per parameter pair, there is a folder which contains one file per gridpoint.
     * Foldertree[0][1] = "OmegamSigma8" would be an example;
     * TrueParamfile is the name of the file which stores whatever theory output is needed to compute the posterior
     * at the true parameter point (which might not be a grid point).
     */
    virtual void AssignTheoryToGridFromFiles(string TreePrefix, vector<vector<string> > Foldertree, vector<string> fileprefixes, string TrueParamFile);
    
    
    virtual void AssignTheoryToGridFromFunction(int p1, int p2, int gp_1, int gp_2);
    
    /* INTERNAL FUNCTIONS. UNLIKELY TO BE MODIFIED AND HENCE NOT VIRTUAL */
    
    
    /**
     * Computes the posterior for the parameter pair p1,p2, and stores it in
     * AllPosteriors. This function is NOT virtual, to avoid the failure mode
     * that users might try to normalze the posterior on the grids, but then
     * forget to divide the posterior at the true parameter point. If that happens
     * then the implemented routine to measure coverages will fail: the coverage
     * is counted by comparing the height of the true point with the height of
     * the confidence contours.
     */
    //Careful: This function will go through OpenMP, do not introduce data races.
    void Posterior_on_Grid(int const p1, int const p2, int const dat_o_sim);
    
    /*
     * Specifies the default contour levels.
     */
    void RequestTheseContourLevels();
    
    double Posterior_at_TrueParams(int const dat_o_sim);
    
    void ComputeCoverage();
    void TwoDimensionalCase(int const p1, int const p2);
    void OneDimensionalCase(int const p1);
    void GridToVector(int const p1, int const p2, int dat_o_sim, vector<double> & output);
    void PosteriorsToFile(int p1, int p2, int dat_o_sim, string  name);
    
    /**
     * This function computes where the confidence contours lie on a posterior.
     * testprob is a vector that stores the posterior likelihood evaluated at
     * as many points as possible. Within the code, it reads in the posterior
     * values from the Posterior_on_Grid(xxx).Confidence contours cut the 
     * posterior at constant height (because they run along isocontours).
     * heights is the vector which stores these values. To be specific, the 
     * function uses vector,double>contourlevels (see above, it is a member variable),
     * and heights are the posterior heights corresponding to these levels.
     * If the string safeto is given, then the contour levels are written to file.
     */
    void ConfidenceContours(vector<double> const testprob, vector<double> & heights, string saveto="");
    
    void pythoncontours( int xcount, int ycount, string executablename, string datafile, int p1, int p2);
    
    void Checks();

    void FreeMemory();

    /**
     * i selects the ith grid point. p the pth param
     * on which delta (the grid spacing) depends.
     */
    double CoordToPhysical(int const p, int const i);
    
    
    /*Member functions for output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    
    /**
     * After the coverage has been measured by
     * analyzing simulations, the results are
     * used to correct the coverage of the post-
     * erior which uses the real data. This function
     * produces parameter plots with the original
     * data analysis setup. It then provides plots
     * with corrected coverage.
     */
    void ApplyCoverageToDataPosterior();
    
    void CoverageTableToFile(int const p1, int const p2, string filename);
    
    void EnglishSummary(string filename="");
    
    
    /*
     * Destructor.
     */
    ~Cover();

};

#endif
