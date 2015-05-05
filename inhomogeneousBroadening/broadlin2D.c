#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

long numTrials = 1000000; /*The size of the ensemble of emitters*/
double R0 = 1.0; /*The minimal defect<-->emitter distance*/
double R = 2000.0; /* The size of the crystal */
long N = 16; /*The number of defects in the crystal*/
//double A[6] = {4.59, 4.59, -4.31, 10.67, 2.63, -2.63}; /*Coupling matrix*/
double A[6] = {0.001, 0.001, 0.001, 1.0, 1.0, 1.0}; /*Coupling matrix*/
/*     ( A[0]    A[3]    A[4] )
A =   ( A[3]    A[1]    A[5] )
        ( A[4]    A[5]    A[2] )
*/
double defStrength = 100000.0; /*Defect's strength*/
double *shift = NULL; /*Array of shifts. Each element corresponds to a single emitter*/

double x = 0.0, y = 0.0, z = 0.0; /*Random coordinates of a defect*/
double r = 0.0;
double strain[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}; /*Strain tensor*/


int main (int argc, char *argv[])
{

  long i = 0;
  long j = 0;

  FILE *shiftsOut = NULL; /*Save the calculated shift into this file*/

  char shiftFileName[2048] = "./";

  shift = (double *) calloc(numTrials, sizeof(double));
  if (shift==NULL)
    {
      printf("Unable to allocate memory for the shift array.\n");
      exit(1);
    }

  if (argc != 12)
    {
      printf("The program requires 11 parameters. Only %d given. \n Falling back to default values.\n", argc-1);
    }
  else
    {
      numTrials = atol(argv[1]);
      R0 = atof(argv[2]);
      R = atof(argv[3]);
      N = atol(argv[4]);
      A[0] = atof(argv[5]);
      A[1] = atof(argv[6]);
      A[2] = atof(argv[7]);
      A[3] = atof(argv[8]);
      A[4] = atof(argv[9]);
      A[5] = atof(argv[10]);
      defStrength = atof(argv[11]);
    }
  strcat(shiftFileName, "pointDefects_2D");
  for (i = 1; i < argc; i++)
    {
      strcat(shiftFileName, "_");
      strcat(shiftFileName, argv[i]);
    }
  strcat(shiftFileName, "\0");
  /*Initializing the random number generator*/
  srand(time(0));

  for (i = 0; i < numTrials; i++)
    {
      memset(strain, 0, 9*sizeof(double));/*Fast clear the stain tensor*/
      for (j=1; j <= N; j++)
        {
          x = 2*R*( (double) rand()/RAND_MAX - 0.5 );
          y = 2*R*( (double) rand()/RAND_MAX - 0.5 );
          r = sqrt(x*x+y*y);
          if ( (r >= R0) && (r <=R) )
            {
              strain[0][0] = strain[0][0] + 1/pow(r,3)*(1-3*pow(x, 2)/pow(r, 2));
              strain[0][1] = strain[0][1] + (-3*x*y/pow(r, 5));
              strain[0][2] = 0.0;

              strain[1][0] = strain[0][1];
              strain[1][1] = strain[1][1] + 1/pow(r,3)*(1-3*pow(y, 2)/pow(r, 2));
              strain[1][2] = 0.0;

              strain[2][0] = strain[0][2];
              strain[2][1] = strain[1][2];
              strain[2][2] = strain[2][2] + 1/pow(r,3);
              shift[i] = defStrength*(A[0]*strain[0][0]+A[1]*strain[1][1]+A[2]*strain[2][2]+2*A[3]*strain[0][1]+2*A[4]*strain[0][2]+2*A[5]*strain[1][2]);
            }
        }    
    }
  /*Save data to file*/
  shiftsOut = fopen(shiftFileName, "w");
  if (shiftsOut==NULL)
    {
      printf("Can not open file. Exiting\n");
      exit(1);
    }
  for (i=0; i < numTrials; i++)
    {
        fprintf(shiftsOut, "%1.20f\n ", shift[i]);
    }
  fclose(shiftsOut); 
  free(shift);
  return 0;
}
