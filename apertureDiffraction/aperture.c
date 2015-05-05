#include <stdlib.h>
#include <math.h>
#include <sys/time.h> // for enhaced randomness
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

#define r1 3.8317 // first zero of Bessel function J_1
#define norm1 2.6318 // normalization for r1

#define r2 7.0156 // second zero of Bessel function J_1
#define norm2 2.8686

#define r3 10.1735//3rd zero
#define norm3 2.9457

#define rX r1 // shortcut to quickly switch between r1, norm1 etc.
#define normX norm1

#define r 0.50//radius of the diffraction spot as ratio r/R
#define MAP_SIZE 400

double R = 0.8;
double w = 2.0; // laser spot e-radius
double g = 60; //saturation parameter;
double norm = 0; // normalization factor for the signal from a single emitter
/* ---------- Monte Carlo Integration Variables------------------*/
size_t dimension = 2;//dimension of the integral in Monte-Carlo
size_t calls = 5000;
double xlower[2] = {0, 0};
double xupper[2] = {2*M_PI, 1.1*r};
gsl_rng *randomizer_MC;
gsl_monte_miser_state *miser_state;

/*----Point Spread Functions. Intensity Distribution Due To Diffraction---*/
double I_Bessel(double rho);
double I_Flattop(double rho);
double I_Partyhat(double rho);
double I_Parabola(double rho);
double I_ErrFuncFamily(double rho);

double Excitation(double *xy);// Excitation of the sample. 
double fraction(double *xy);
double f(double *x, size_t dim, void *params);

int
main (void)
{

  const gsl_rng_type * T;
  gsl_rng * randomizer; 

  double signal = 0; // signal from the current point on the map
  double S1 = 0;

  unsigned int i, n;
  double xy[2]; // random position of the emitting center

  double S_diffr = 0;
  double S_geom = 0;
  double delta = 0;
  struct timeval tv;
  long seed;

  size_t Map_Size = MAP_SIZE;
  double SGeom_Set[MAP_SIZE];
  double SDiffr_Set[MAP_SIZE];
  double Delta_Set[MAP_SIZE];
  double N_Avg = 5000;
  double N_emit; // number of emitting centers in the working area
  double N_inside = 0;
  // statistical characteristics 
  double SGeom_Mean, SGeom_Sigma;
  double SDiffr_Mean, SDiffr_Sigma;
  double Delta_Mean, Delta_Sigma;
  double BetaGeom;
  double BetaDiffr;

  miser_state = gsl_monte_miser_alloc(dimension);

  {//setting up randomizer
    gsl_rng_env_setup();
    T = gsl_rng_ranlux389;
    randomizer = gsl_rng_alloc (T);// for generating centers' coordinates
    randomizer_MC= gsl_rng_alloc (T); // for Monte-Carlo integration
    gettimeofday(&tv, NULL);
    seed = tv.tv_usec;
    gsl_rng_set(randomizer, seed);
    gettimeofday(&tv, NULL);
    seed = tv.tv_usec;
    gsl_rng_set(randomizer_MC, seed);
  }
  // to warm-up randomizer
  xy[0] = gsl_ran_flat(randomizer, -(R+r), +(R+r));
  xy[1] = gsl_ran_flat(randomizer, -(R+r), +(R+r));
  /* calculating total signal from a single center in the center of
     the aperture. We need to increase the size of the aperture for a moment not to
     loose the signal if r>R */
  xy[0] = 0.0;
  xy[1] = 0.0;
  R = 2*r;
  S1 = fraction(xy);
  R = 1.0000; // closing the aperture to it's initial size
  printf("Total signal: %5.4f\n", S1);

  for (n = 0; n < Map_Size; n++)
    {
      S_diffr = 0;
      S_geom = 0;
      delta = 0;      
      N_emit = gsl_ran_poisson(randomizer, N_Avg);
      // calculating the signals at a given point on the map
      for (i = 1; i <= N_emit; i++)
	{
	  xy[0] = gsl_ran_flat(randomizer, -(R+r), +(R+r));
	  xy[1] = gsl_ran_flat(randomizer, -(R+r), +(R+r));	  
	  signal = S1*Excitation(xy);// total signal from the center at this position {xi, yi}
      
	  if ( (fabs(xy[0]) < R) && (fabs(xy[1]) < R) )
	    // the emitter is inside the aperture
	    {
	      S_geom = S_geom + signal;
	      N_inside = N_inside + 1;
	      if ( (fabs(xy[0]) < (R-r) ) && (fabs(xy[1]) < (R-r) ) )
		// the center is inside the "core" -- no need to integrate
		{
		  S_diffr = S_diffr+signal;
		}
	      else
		{
		  S_diffr = S_diffr + signal*fraction(xy);
		}
	    }
	  else
	    {
	      S_diffr = S_diffr + signal*fraction(xy);
	    }
	}
      delta = S_diffr/S_geom;
      SGeom_Set[n] = S_geom;
      SDiffr_Set[n] = S_diffr;
      Delta_Set[n] = delta;
    }

  {
    gsl_rng_free(randomizer);
    gsl_rng_free(randomizer_MC);
    gsl_monte_miser_free(miser_state);
  }

  {// Extracting the statistical parameters
    SGeom_Mean = gsl_stats_mean(SGeom_Set, 1, MAP_SIZE);
    SGeom_Sigma = gsl_stats_sd(SGeom_Set, 1, MAP_SIZE);

    SDiffr_Mean = gsl_stats_mean(SDiffr_Set, 1, MAP_SIZE);
    SDiffr_Sigma = gsl_stats_sd(SDiffr_Set, 1, MAP_SIZE);

    Delta_Mean = gsl_stats_mean(Delta_Set, 1, MAP_SIZE);
    Delta_Sigma = gsl_stats_sd(Delta_Set, 1, MAP_SIZE);

    BetaGeom  = ( N_Avg/( (1+r)*(1+r) ) )*SGeom_Sigma*SGeom_Sigma/(SGeom_Mean*SGeom_Mean);
    BetaDiffr = N_Avg*SDiffr_Sigma*SDiffr_Sigma/(SDiffr_Mean*SDiffr_Mean);
  }

  {// Displaying the results
    printf("Beta Geometrical: %5.4f\n", BetaGeom);
    printf("Beta Diffractional: %5.4f\n", BetaDiffr);
    printf("Average Delta: %5.4f\n", Delta_Mean);
    printf("Delta STD: %5.4f\n", Delta_Sigma);
  }

  return 0;
}

double I_Bessel(double rho)
// Intensity distribution in the Point Spread Function
{
  double result = 0;
  double x = rho*rX/r;
  if (rho < r)
    {
      if (x == 0)
	{
	  result = (rX*rX)/(4*normX*r*r);
	}
      else
	{
	  result = (rX*rX)/(normX*r*r)*(gsl_sf_bessel_J1(x)/x)*(gsl_sf_bessel_J1(x)/x);
	}
    }
  else
    result = 0;
  return result;
}

double I_Flattop(double rho)
//Normalized flat-top circle function
{
  double result = 0;
  if (rho < r)
    result = 1/(M_PI*r*r);
  else
    result = 0;
  return result;
}

double I_Partyhat(double rho)
{
  double result = 0;
  if (rho < r)
    result = (1-rho/r)/0.2619;//normalized
  else
    result = 0;
  return result;
}

double I_Parabola(double rho)
{
  double result = 0;
  if (rho < r)
    result = (1-(rho*rho)/(r*r))/0.3927;//normalized
  else
    result = 0;
  return result;
}

double I_ErrFuncFamily(double rho)
{
  double result = 0;
  double scale = 5; // the bigger the scale --- the closer the the flat-top
  if (rho < r)
    result = erf(scale*(r-rho))/0.4981;
  else
    result = 0;
  return result;
}  

double Excitation(double *xy)
/* Distribution of the intensities of emitting centers due to
   the excitation profile and saturation effects*/
{
  return exp( (xy[0]*xy[0]+xy[1]*xy[1])/(w*w) )/(1+g*exp( (xy[0]*xy[0]+xy[1]*xy[1])/(w*w) ));
}

double f(double *x, size_t dim, void *params)
// the function being integrated in polar coordinates I(r)*r
{
  double result = 0;
  double xj, yj;
  double *xiyi = (double *) params;
  //calculating the coordinates with respect to the aperture
  xj = xiyi[0] + x[1]*cos(x[0]);
  yj = xiyi[1] + x[1]*sin(x[0]);
  
  if ( ( fabs(xj) < R ) && ( fabs(yj) < R ) )
    result = I_Bessel(x[1])*x[1];
  else
    result = 0;
  return result;
}

double fraction(double *xy)
{
  gsl_monte_function F = {&f, dimension, (void *)xy};
  double result, error;

  gsl_monte_miser_init(miser_state);
  gsl_monte_miser_integrate(&F, xlower, xupper, dimension,
			    calls, randomizer_MC, miser_state,
			    &result, &error);
  //printf("%5.4f\n", result);
  return result;
}
