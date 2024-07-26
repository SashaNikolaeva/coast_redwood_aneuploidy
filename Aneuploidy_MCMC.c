// This code is written by Rasmus Nielsen
// Feel free to use as you like, but it is always nice if you cite our papers

// uniform Dirichlet prior for alpha
// truncated geometric prior for k, 0<=k<=12, with Pr(k=6) = p*c, Pr(k=5)=Pr(k=7) = c*p(1-p)/2, etc. c = \sum_{i=0}^6 p(1-p)^i, 
// U[0,1] hyperprior for p.
// multinomial likelihood as a product over all i and j of (a_i*k_{ij}/(\sum_{v=1}^k a_v*k_{vj}))^n_{ij}
// reflected exponential to marginally update each alpha_i
// reflected exponential to marginally update p
// random walk to update k_{ij} marginally

#include <stdlib.h>
#include <stdio.h>
#include <math.h> 

FILE *infile, *parfile, *outfile, *trace_alpha_file, *trace_p_file;
#define VERBOSE 1
#define EXPECTED_PLOIDY 6 // must be even

int log_of_one; // constant we are going to need and that we start by calculating in Main

// Random number generator
/*U(0,1): AS 183: Appl. Stat. 31:188-190
Wichmann BA & Hill ID.  1982.  An efficient and portable pseudo-random number generator.  Appl. Stat. 31:188-190 x, y, z are any numbers in the range 1-30000.  Integer operation up to 30323 required.

Suggested to me by Z. Yang who also provided me with the source code used here. */

int static z_rndu = 137;

void SetSeed(int seed) {
  z_rndu = 170 * (seed % 178) + 137;
}

double uniform() {
  static int x_rndu = 11, y_rndu = 23;

  double r;
  x_rndu = 171 * (x_rndu % 177) -  2 * (x_rndu / 177);
  y_rndu = 172 * (y_rndu % 176) - 35 * (y_rndu / 176);
  z_rndu = 170 * (z_rndu % 178) - 63 * (z_rndu / 178);
  if (x_rndu < 0) x_rndu += 30269;
  if (y_rndu < 0) y_rndu += 30307;
  if (z_rndu < 0) z_rndu += 30323;
  r = x_rndu / 30269.0 + y_rndu / 30307.0 + z_rndu / 30323.0;
  return (r - (int)r);
}

// simulates a proposed update using exponential RV constrained to be in [0,max] and reflected at boundaries
double exponential_reflected(double rate, double c, int max)
{
  double a;

  a = uniform();
  if (a > 0.5) c = c - log(2.0 * (a - 0.5)) / rate;
  else c = c + log(2.0 * a) / rate;
  while (c < 0 || c > max)
  {
    if (c < 0) c = -c;
    if (c > max) c = 2.0 * max - c;
  }
  return c;
}

// This function calculates the log likelihood given the latent values of k and the alpha vector
double loglike(double *alpha, int numchrom, int numindivid, int **ploidy_matrix, int **data, double fudgefactor)
{
  int i, j;
  double p, sum;

  p = 0;
  for (j = 0; j < numindivid; j++) {
    sum = 0;
    for (i = 0; i < numchrom; i++)
      sum = sum + alpha[i] * ploidy_matrix[i][j];
    for (i = 0; i < numchrom; i++) {
      p = p + (double)data[i][j] * log(alpha[i] * (double)ploidy_matrix[i][j] / sum) / fudgefactor;
    }
  }
  return p;
}

// This function calculates the partial log likelihood, only for individual 'ind', given the latent values of k and the alpha vector
double ind_loglike(double *alpha, int numchrom, int ind, int numindivid, int **ploidy_matrix, int **data, double fudgefactor)
{
  int i;
  double p, sum;

  if (ind < 0 || ind >= numindivid) { printf("error in ind_loglike"); exit(-1); }
  p = 0;
  sum = 0;
  for (i = 0; i < numchrom; i++)
    sum = sum + alpha[i] * ploidy_matrix[i][ind];
  for (i = 0; i < numchrom; i++)
    p = p + (double)data[i][ind] * log(alpha[i] * ploidy_matrix[i][ind] / sum) / fudgefactor;
  return p;
}

// this function updates the alpha vector
// alphanew is scrap memory for temporarily storing the new value of alpha
int update_alpha(double *alpha, double *alphanew, double rate, int numchrom, int numindivid, int chrom, int **ploidy_matrix, int **data, double fudgefactor)
{
  double a, r, diffac; 
  int i;

  a = exponential_reflected(rate, alpha[chrom], 1); // new proposed value of alpha for chromosome chrom
  diffac = (a - 1.0) / (alpha[chrom] - 1.0); // this is the multiplicative factor that all the other entries in alpha must be updated with for the alphas to sum to one
  for (i = 0; i < numchrom; i++) { // now make the new alpha vector
    if (i == chrom) alphanew[i] = a;
    else alphanew[i] = alpha[i] * diffac;
  }
  // calculate acceptance probability. Because of symmetric proposal probability and uniform dirichlet, prior or proposal prob does not enter
  r = loglike(alphanew, numchrom, numindivid, ploidy_matrix, data, fudgefactor) - loglike(alpha, numchrom, numindivid, ploidy_matrix, data, fudgefactor);
  if (r > log_of_one || log(uniform()) < r) return 1; // if accepted - switch pointers between alpha and alphanew
  else return 0;
}

double log_prior_ratio_k(int oldk, int newk, double p)
{
  if (oldk == EXPECTED_PLOIDY) return  log((1.0 - p) / 2.0); 
  else if (newk == EXPECTED_PLOIDY) return log(2.0 / ((1.0 - p)));
  else if (abs(oldk - EXPECTED_PLOIDY) < abs(newk - EXPECTED_PLOIDY)) return log(1.0 - p);
  else return log(1.0 / (1 - p)); 
}

// this function updates k for individual ind, chromosome chrom
void update_k(double *alpha, int ind, double p, int chrom, int numchrom, int numindivid, int maxploidy, int **ploidy_matrix, int **data, double fudgefactor)
{
  double oldk, r;

  r = -ind_loglike(alpha, numchrom, ind, numindivid, ploidy_matrix, data, fudgefactor); // likelihood component of acceptance probability for the old state.
  oldk = ploidy_matrix[chrom][ind];
  if (uniform() < 0.5) ploidy_matrix[chrom][ind] = oldk + 1; // we propose an update up or down with equal probability
  else ploidy_matrix[chrom][ind] = oldk - 1;  
  if (ploidy_matrix[chrom][ind] == -1) ploidy_matrix[chrom][ind] = maxploidy; // We arrange it on a ring to deal with boundaries
  else if (ploidy_matrix[chrom][ind] > maxploidy)  ploidy_matrix[chrom][ind] = 0; 
  if (ploidy_matrix[chrom][ind] == 0 && data[chrom][ind] > 0) ploidy_matrix[chrom][ind] = oldk; // done to avoid some underflow problems when ploidy is zero but there are reads
  else {
    r = r + ind_loglike(alpha, numchrom, ind, numindivid, ploidy_matrix, data, fudgefactor); // likelihood component of acceptance probability for the new state. No proposal component as it is symmetric
    r = r + log_prior_ratio_k(oldk, ploidy_matrix[chrom][ind], p); // prior part of update probability
    if (r < log_of_one  && log(uniform()) > r) // if not accepted - switch back
      ploidy_matrix[chrom][ind] = oldk;
  }
}

double update_p(double p, double rate, int numchrom, int numindivid, int **ploidy_matrix)
{
  int i, j, countdif = 0;
  double pnew, r;

  pnew = exponential_reflected(rate, p, 1);
  for (j = 0; j < numindivid; j++)
    for (i = 0; i < numchrom; i++)
      countdif += abs(EXPECTED_PLOIDY - ploidy_matrix[i][j]);
  r = (numchrom * numindivid) * log(pnew / p) + countdif * log((1.0 - pnew) / (1.0 - p)); // contribution to augmented likelihood, proposal symmetric and prior is uniform
  if (r > log_of_one || log(uniform()) < r) // if  accepted  - again all is done in log space to avoid underflow/overflow
    return pnew;
  else return p;
}

void mcmc(int **data, int numchrom, int numindivid, int rep, int burnin, double rate_update_p, double rate_update_alpha, char **names, char **chromosomenames)
{
  int i, j, r, k, **ploidy_matrix, maxploidy;
  double *alpha, *alphanew, *store, ***ploidyprob, *avealpha, p = 0.001, fudgefactor;

  maxploidy = 2 * EXPECTED_PLOIDY; 
  // First allocate all kind of memory needed
  alpha = malloc(numchrom * sizeof(double));
  alphanew = malloc(numchrom * sizeof(double));
  avealpha = malloc(numchrom * sizeof(double));
  for (i = 0; i < numchrom; i++)
    alpha[i] = 1.0 / (double)(numchrom);
  ploidy_matrix = malloc(numchrom * sizeof(int *));
  for (i = 0; i < numchrom; i++)
    ploidy_matrix[i] = malloc(numindivid * sizeof(int));
  for (i = 0; i < numchrom; i++)
    for (j = 0; j < numindivid; j++)
      ploidy_matrix[i][j] = EXPECTED_PLOIDY;
  ploidyprob = malloc(numchrom * sizeof(double **));
  for (i = 0; i < numchrom; i++) {
    avealpha[i] = 0;
    ploidyprob[i] = malloc(numindivid * sizeof(double *));
    for (j = 0; j < numindivid; j++) {
      ploidyprob[i][j] = malloc((maxploidy + 1) * sizeof(double));
      for (k = 0; k <= maxploidy; k++)
        ploidyprob[i][j][k] = 0.0;
    }
  }
  fudgefactor = 1.0; // set to 1 for regular mcmc

  // Doing some updates of alpha first
  for (r = 0; r < rep / 20; r++) {
    for (i = 0; i < numchrom; i++) {
      if(update_alpha(alpha, alphanew, rate_update_alpha, numchrom, numindivid, i, ploidy_matrix, data, fudgefactor)) {
        store = alpha; // swap pointers if update is accepted
        alpha = alphanew;
        alphanew = store;
      }
    }
  }

  // Now start the 'rep' number of MCMC iterations
  if (VERBOSE)  printf("Iteration number:\n");
  for (r = 0; r < rep; r++) {
    if (VERBOSE) { if (1000 * (r / 1000) == r) { printf("%i\n", r); }}
    for (i = 0; i < numchrom; i++) {
      for (j = 0; j < numindivid; j++)
        update_k(alpha, j, p, i, numchrom, numindivid, maxploidy, ploidy_matrix, data, fudgefactor); // updates the ploidy
      if(update_alpha(alpha, alphanew, rate_update_alpha, numchrom, numindivid, i, ploidy_matrix, data, fudgefactor)) { // updates alpha
        store = alpha; // swap pointers if update is accepted
        alpha = alphanew;
        alphanew = store;
      }
    }
    p = update_p(p, rate_update_p, numchrom, numindivid, ploidy_matrix); // updates p

    // Record trace data
    if (r >= burnin) {
      for (i = 0; i < numchrom; i++) {
        fprintf(trace_alpha_file, "%.8f ", alpha[i]);
      }
      fprintf(trace_alpha_file, "\n");
	  fprintf(trace_p_file, "%.8f\n", p);
      
    }
  }

  // Print the results
  for (j = 0; j < numindivid; j++) {
    fprintf(outfile, "Individual %s:\n", names[j]);
    for (k = 0; k <= maxploidy; k++)
      fprintf(outfile, "\t%i", k);
    fprintf(outfile, "\n");
    for (i = 0; i < numchrom; i++) {
      fprintf(outfile, "%s\t", chromosomenames[i]);
      for (k = 0; k <= maxploidy; k++)
        fprintf(outfile, "%.3lf\t", ploidyprob[i][j][k] / (double)(rep - burnin));
      fprintf(outfile, "\n");
    }
  }
  fprintf(outfile, "\nExpected dosages per chromosome:\n");
  for (i = 0; i < numchrom; i++)
    fprintf(outfile, "%s %.4lf\n", chromosomenames[i], avealpha[i] / (double)(rep - burnin));
  
  // Free the memory allocated
  for (i = 0; i < numchrom; i++)
    free(ploidy_matrix[i]);
  free(ploidy_matrix);
  for (i = 0; i < numchrom; i++) {
    for (j = 0; j < numindivid; j++)
      free(ploidyprob[i][j]);
    free(ploidyprob[i]);
  }
  free(ploidyprob);
  free(alpha);
  free(alphanew);
  free(avealpha);
}

int main(int argc, char *argv[])
{
  char inname[30], outname[30], parname[30], **names, **chromosomenames;
  int i, j, **data, numreps, burnin, seed, numchrom, numindivid;
  double rate_update_p, rate_update_alpha;

  log_of_one = log(1.0);

  if (argc < 4) { printf("Specify name of the infile with data (first), outfile (second) and parameterfile (third)"); exit(-1); }
  sprintf(inname, "%s", argv[1]);
  sprintf(outname, "%s", argv[2]);
  sprintf(parname, "%s", argv[3]);

  // This code reads in parameters
  if (NULL == (parfile = fopen(parname, "r"))) {
    puts("Cannot open parameterfile!");
    exit(-1);
  }
  fscanf(parfile, "%i %i %i %lf %lf", &numreps, &burnin, &seed, &rate_update_p, &rate_update_alpha);
  if (burnin >= numreps) { printf("The number of iterations must be larger than the burnin\n"); exit(-1); }
  fclose(parfile);

  // This code reads the data
  if (NULL == (infile = fopen(inname, "r"))) {
    puts("Cannot open infile!");
    exit(-1);
  }
  fscanf(infile, "%i %i\n", &numchrom, &numindivid);
  data = malloc(numchrom * sizeof(int *));
  for (i = 0; i < numchrom; i++)
    data[i] = malloc(numindivid * sizeof(int));
  chromosomenames = malloc(numchrom * sizeof(char *));
  for (i = 0; i < numchrom; i++)
    chromosomenames[i] = malloc(100 * sizeof(char));
  for (i = 0; i < numchrom; i++)
    fscanf(infile, "%s", chromosomenames[i]);
  names = malloc(numindivid * sizeof(char *));
  for (j = 0; j < numindivid; j++)
    names[j] = malloc(100 * sizeof(char));
  for (j = 0; j < numindivid; j++) {
    fscanf(infile, "%s", names[j]);   
    if (VERBOSE) printf("Found data for %s\n", names[j]);
    for (i = 0; i < numchrom; i++)
      fscanf(infile, "%i", &data[i][j]);
  }
  fclose(infile);
  
  // Opens the outfile so it is ready for input
  if (NULL == (outfile = fopen(outname, "w"))) {
    puts("Cannot open outfile!");
    exit(-1);
  }

  // Open trace files for alpha and p
  trace_alpha_file = fopen("/space/s1/sashanikolaeva/Scripts/Rasmus_MCMC/trace_alpha.txt", "w");
  if (trace_alpha_file == NULL) {
    printf("Cannot open trace file for alpha values.\n");
    exit(-1);
  }
  trace_p_file = fopen("/space/s1/sashanikolaeva/Scripts/Rasmus_MCMC/trace_p.txt", "w");
  if (trace_p_file == NULL) {
    printf("Cannot open trace file for p values.\n");
    exit(-1);
  }

  // MCMC starts here
  SetSeed(seed);
  mcmc(data, numchrom, numindivid, numreps, burnin, rate_update_p, rate_update_alpha, names, chromosomenames);

  // Cleaning up - need to also do chromosome names and sample names
  fclose(outfile);
  fclose(trace_alpha_file);
  fclose(trace_p_file);
  
  // Free the memory allocated
  for (i = 0; i < numchrom; i++)
    free(data[i]);
  free(data);
}
