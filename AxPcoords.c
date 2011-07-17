/*  Axelerated Pcoords Implementation                                           */
/*  Copyright  January 2007 by Alexandros Stamatakis                            */
/*  All rights reserved.                                                        */
/*                                                                              */
/*  This is a highly optimized and parallelized C porting of the original       */
/*  Pcoords Fortran code by Pierre Legendre available at                        */
/*                                                                              */
/*  http://www.bio.umontreal.ca/Casgrain/en/labo/parafit.html                   */
/*                                                                              */
/*  and published as:                                                           */
/*                                                                              */ 
/*  Legendre, P., Y. Desdevises and E. Bazin. 2002.                             */ 
/*  A statistical test for host-parasite coevolution.                           */
/*  Systematic Biology 51(2): 217-234.                                          */
/*                                                                              */
/*                                                                              */
/*  This code may be used and modified for non-commercial purposes              */
/*  but redistribution in any form requires written permission.                 */
/*  Please contact:                                                             */
/*                                                                              */
/*  Alexandros Stamatakis                                                       */       
/*  Swiss Federal Institute of Technology                                       */
/*  School of Computer & Communication Sciences                                 */
/*  Laboratory for Computational Biology and Bioinformatics (LCBB)              */
/*  STATION 14                                                                  */
/*  CH-1015 Lausanne, Switzerland                                               */
/*                                                                              */
/*  Tel:   +41 21 69 31392 (Office)                                             */
/*         +41 22 54 80003 (SkypeIn)                                            */
/*         +41 796115849   (Mobile)                                             */
/*  Skype: stamatak                                                             */
/*  Email: Alexandros.Stamatakis@epfl.ch                                        */
/*  WWW:   icwww.epfl.ch/~stamatak                                              */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#ifdef _USE_GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#endif
  
#ifdef _USE_ACML
#include <acml.h>
#endif

#ifdef _USE_MKL
#include <mkl.h>
#endif

#ifdef _WIN32
#include <time.h>
#endif



#define EIGEN_EPSILON 0.000005

/* The constant above could be adapted appropriately */

double gettime()
{
  #ifndef _WIN32
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
  #else
  clock_t time=clock();
  return ((double)time) / ((double)CLOCKS_PER_SEC);
  #endif
}




/* Interface Stuff ***************************************************************************************************************************/

static int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"r");
  
  if(fp) 
    {
      res = 1;
      fclose(fp);
    }
  else 
    res = 0;
       
  return res;
} 

typedef struct {
  int n;
  int transposed;
  char inputFileName[1024];  
  char logFileName[1024];
  char outFileName[1024]; 
} parameters;



static void printHelp()
{
  printf("\n\n\n\n");
  
  printf("AxPcoords version 1.0, February 2007 by Alexandros Stamatakis\n");
  printf("Please also consult the Manual\n");
  printf("To report bugs send an email to Alexandros.Stamatakis@epfl.ch\n\n\n");
  
  printf("AxPcoords[GSL|ACML|MKL] -f inputMatrixName -n rowNumber [-h] [-t]\n");  
  printf("\n");
  printf("       -f     Specify file Name of quadratic input matrix\n");  
  printf("\n");
  printf("       -n     Specify number of Rows in quadratic input Matrix\n");
  printf("\n");
  printf("       -h     Display this help message\n");  
  printf("\n");
  printf("       -t     Specify if the output matrix shall be transposed\n");
  printf("              DEFAULT: OFF\n");   
  printf("\n\n\n\n");
}


static void get_args(int argc, char *argv[], parameters *params)
{
  int i, k, res;
  int number;
#define NUM_OPT 4
  char *options[NUM_OPT]= {"-f", "-n", "-t", "-h"};
  char fileName[1024];
  char runName[1024];
  int set[NUM_OPT];  

  /* init */

  params->n = 0;
  params->transposed = 0;
 

  /***************/

  for(i = 0; i < NUM_OPT; i++)
    set[i] = 0;

  for(i = 1; i < argc; i++)
    {
      int found = 0;
      
      for(k = 0; k < NUM_OPT && !found; k++)	
	if(strcmp(options[k], argv[i]) == 0)	    	    	    
	  found = 1;	      	      	         
      k--;
      
      if(found)
	{	 
	  switch(k)
	    {
	    case 0:
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  strcpy(fileName, argv[++i]);		   
		  
		  if(!filexists(fileName))
		    {
		      printf("file %s does not exist\n", fileName);
		      exit(-1);
		    }
		  strcpy(params->inputFileName, fileName);

		  strcpy(params->logFileName, fileName);
		  strcat(params->logFileName, ".pcoords.log");

		   if(filexists(params->logFileName))
		    {
		      printf("output file %s does already exist\n", params->logFileName);
		      exit(-1);
		    }

		  strcpy(params->outFileName, fileName);
		  strcat(params->outFileName, ".pcoords");
		  
		  if(filexists(params->outFileName))
		    {
		      printf("output file %s does already exist\n", params->outFileName);
		      exit(-1);
		    }

		}
	      strcpy(params->inputFileName, fileName);
	      set[k] = 1;
	      break;
	    case 1:
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  res = sscanf(argv[++i],"%d", &number);
		  if(res == 0 || res == EOF)
		    {
		      printf("argument of %s must be an integer number\n", options[k]);
		      exit(-1);
		    }
		}
	      params->n = number;
	      set[k] = 1;		
	      break;
	    case 2:
	      params->transposed = 1;
	      break;	      	  	      	   
	    case 3:	      
	      printHelp();
	      exit(0);
	      break;	   
	    default:
	      printf("unknown option %s\n", options[k]);
	      exit(-1);
	    }
	}
      else
	{
	  printf("unknown option %s\n", argv[i]);
	  exit(-1);
	}     
    }


  for(i = 0; i < 2; i++)
    {
      if(set[i] == 0)
	{
	  printf("Error option %s must be specified!\n", options[i]);
	  exit(-1);
	}
    }   
}

/************************* Interface end ***************************************************************************************************/

int main (int argc, char *argv[])
{
  unsigned int n, nSquare;
  double *d, *d2;
  double avg, value, nSquareD, nD, trace;
  double *rowAvg, *colAvg;
  FILE *f, *outf, *logf;
  char fileName[1024], outFileName[1024], logFileName[1024];
  double time;
  unsigned int i, j;
  int posEV = 0; 
  parameters *params = (parameters *)malloc(sizeof(parameters));

  time = gettime();  

  get_args(argc, argv, params);

  n = params->n;
 
  strcpy(fileName,    params->inputFileName);
  strcpy(outFileName, params->outFileName);
  strcpy(logFileName, params->logFileName);  
  
  outf = fopen(outFileName, "w");
  logf = fopen(logFileName, "w");

  fprintf(logf, "Distance Matrix Size: %d\n", n);  
  fprintf(logf, "Input file name: %s\n", fileName);

  nSquare = n * n;

  nSquareD = (double)nSquare;
  nD       = 1/((double)n);
   
  d  = (double *)malloc(sizeof(double) * nSquare);
  
  rowAvg = (double *)malloc(sizeof(double) * n);
  colAvg = (double *)malloc(sizeof(double) * n);          

  f = fopen(fileName, "r");  
  for(i = 0; i < nSquare; i++)     
    {
      int v;
      v = fscanf(f, "%lf",&d[i]);    
      if(v == 0)
	{
	  printf("Format Conversion Error while reading Matrix in File %s at position %d \n", fileName, i);
	  exit(-1);
	}
      if(v == EOF)
	{
	  printf("End of File reached while reading Matrix in File %s at position %d\n", fileName, i);
	  exit(-1);
	}
    }
  fclose(f); 

  avg = 0; 
  for(i = 0; i < nSquare; i++)
    {
      value = (-0.5) * d[i] * d[i];
      d[i]  = value;
      avg  += value;
    }  

  avg = avg / nSquareD;
   
  for(i = 0; i < n; i++)
    {
      value = 0;
      for(j = 0; j < n; j++)	
	value += d[i * n + j];	
      rowAvg[i] = value * nD;
    }
  
  for(j = 0; j < n; j++)
    {
      value = 0;
      for(i = 0; i < n; i++)
	{
	  value += d[i * n + j];
	}
      colAvg[j] = value * nD;
    }  
    
  for(i = 0; i < n; i++)
    {       
      for(j = 0; j < n; j++)
	d[i * n + j] =  d[i * n + j]  - colAvg[j] - rowAvg[i] + avg;
    } 

  trace = 0.0;

  for(i = 0; i < n; i++)
    trace += d[i * n + i];

  printf("START TRACE %f\n", trace);
  fprintf(logf, "START TRACE %f\n", trace);
 
  free(rowAvg);
  free(colAvg);

#ifdef _USE_LAPACK
  {
    int info, lwork, posLimit; 
    double *w = (double *)calloc(n, sizeof(double));
    
    #ifdef _USE_ACML
    dsyev('V','U', n, d, n, w, &info);
    #elif _USE_MKL
    double workSize=0.0; // size of work space
    double* work; // work space
    char jobz='V';
    char uplo='U';

    lwork=-1; // i.e., return estimated work space size
    dsyev(&jobz,&uplo, &n, d, &n, w, &workSize, &lwork, &info);
    
    // allocate work space    
    lwork=(int)workSize;
    work=(double*) malloc(sizeof(double) * lwork);
    
    // do the real work
    dsyev(&jobz,&uplo, &n, d, &n, w, work, &lwork, &info);
    
    free(work);
    #endif

    trace = 0.0;
      
    for(i = 0; i < n; i++)
      trace += w[i];
      
    printf("TRACE %f\n", trace);
    fprintf(logf, "EV-TRACE %f\n\n", trace);
   
    fprintf(logf, "Eigenvalues:\n");

    for(j = 1; j <= n; j++)
      fprintf(logf, "%1.5f \n", w[n - j]);
      
    for(j = 0; j < n; j++)
      {       
	double temp;
	if(w[j] > EIGEN_EPSILON)	    
	  posEV++;	    
	
	if(w[j] > 0.0)	  
	  temp = sqrt(w[j]);
	else
	  temp = sqrt(-w[j]);	

	for(i = 0; i < n; i++)
	  d[j * n + i] *= temp;
      }
    
    printf("\n Number of positive Eigenvalues: %d\n", posEV);    
    fprintf(logf, "\n Number of positive Eigenvalues: %d\n", posEV);  

    posLimit = n - posEV - 1;


    if(params->transposed == 1)
      {
	for(i = n - 1; i > posLimit; i--)
	  for(j = 0; j < n; j++)	  	    	   	      
	    fprintf(outf, "%1.5f \n", d[i * n + j]);			      	  
      }
    else
      {
	for(j = 0; j < n; j++)
	  {
	    i = n - 1;
	    
	    while(i > posLimit)	 
	      {
		fprintf(outf, "%1.5f \n", d[i * n + j]);
		i--;
	      }	
	  }	        
      }
  }
#endif

#ifdef _USE_GSL   
  { 
    double sum = 0;
    double *eigenValues, *eigenVectors;
      
    gsl_matrix_view m = gsl_matrix_view_array(d, n, n);
    gsl_vector *eval = gsl_vector_alloc(n);
    gsl_matrix *evec = gsl_matrix_alloc(n, n);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);       
    gsl_eigen_symmv (&m.matrix, eval, evec, w);     
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);    
    eigenValues  = gsl_vector_ptr(eval, 0);
    eigenVectors = gsl_matrix_ptr(evec, 0, 0);
     
    trace = 0.0;

    for(i = 0; i < n; i++)
      trace += eigenValues[i];           

    printf("TRACE %f\n", trace);
    fprintf(logf, "EV-TRACE %f\n\n", trace);
   
    fprintf(logf, "Eigenvalues:\n");

    for(j = 0; j < n; j++)
      fprintf(logf, "%1.5f \n", eigenValues[j]);

    for(j = 0; j < n; j++)
      {
	double temp;

	if(eigenValues[j] > EIGEN_EPSILON)	    
	  posEV++;	    
	
	if(eigenValues[j] > 0.0)	  
	  temp = sqrt(eigenValues[j]);
	else
	  temp = sqrt(-eigenValues[j]);

	
	for(i = 0; i < n; i++)
	  eigenVectors[i * n + j] *= temp;	 	
      }
    
    printf("Number of positive Eigenvalues: %d\n", posEV);             
    fprintf(logf, "\n Number of positive Eigenvalues: %d\n", posEV);


    if(params->transposed == 1)
      {
	for(i = 0; i < posEV; i++)	  
	  for(j = 0; j < n; j++)	        
	    fprintf(outf, "%1.5f \n", eigenVectors[j * n + i]);         	  
      }
    else
      {
	for(j = 0; j < n; j++)	  
	    for(i = 0; i < posEV; i++)	        
	      fprintf(outf, "%1.5f \n", eigenVectors[j * n + i]);         	   
      }
  }
#endif
    
  printf("Results have been written to file: %s\n", outFileName);
  fprintf(logf, "Results have been written to file: %s\n", outFileName);


  printf("\nEXECUTION TIME %f\n", gettime() - time);
  fprintf(logf, "\nEXECUTION TIME %f\n", gettime() - time);

  fclose(logf);
  fclose(outf);
  
 
  return 0;
}
