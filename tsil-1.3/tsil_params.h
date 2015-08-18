#ifndef TSIL_PARAMS_H
#define TSIL_PARAMS_H

/* This file contains parameters associated with the numerical
   integration in TSIL */

/* The macros TSIL_PRECISION_GOAL, TSIL_NSTEPS_START,
   TSIL_NSTEPS_MAX_CON, TSIL_NSTEPS_MAX_VAR, and TSIL_NSTEPS_MIN are
   used by TSIL_SetParameters to set the corresponding members of the
   data struct. To reset them at run time, use the function
   TSIL_ResetStepSizeParams after calling TSIL_SetParameters */

#define TSIL_PRECISION_GOAL 1e-12
/* For each RK step, require that the estimated error for each
   dependent variable is less than TSIL_PRECISION_GOAL times the
   increment of that dependent variable for that step. Otherwise, the
   step size is reduced and the step is retried.  This requirement is
   ignored if the step size gets too small, or if the estimated error
   is less than TSIL_EPSILON [DGR - TSIL_TOL?] times the absolute value of the dependent
   variable. If TSIL_PRECISION_GOAL is set to 0, the step sizes will
   be determined by TSIL_NSTEPS_MAX_CON and TSIL_NSTEPS_MAX_VAR below.
   If TSIL_PRECISION_GOAL is set to a large number, the step sizes
   will tend towards that given by TSIL_NSTEPS_MIN below. */

#define TSIL_NSTEPS_START 500
/* For each leg of the contour of the RK integration, the initial step
   size is chosen so that there would be 500 total steps if the step
   size did not change. */

#define TSIL_NSTEPS_MAX_CON 10000
#define TSIL_NSTEPS_MAX_VAR 10000
#define TSIL_NSTEPS_MIN     100
/* The minimum step size for a leg of the contour with dimensionless
   rescaled independent variable length L is given by:

   L/(N_STEPS_MAX_CON + L * N_STEPS_MAX_VAR).

   The maximum step size is L/NSTEPS_MIN. */

#define IM_DISPL 0.20L
/* Contour displacement in Im(s) direction to avoid thresholds, in
   terms of dimensionless rescaled variable. */

#define THRESH_CUTOFF 0.025L
/* Maximum distance considered "close" to a threshold, in terms of
   dimensionless rescaled variable. */

/* DGR modified from TSIL_EPSILON */
#define SINIT 10.0*TSIL_TOL
/* Starting point for integration when s=0 is (or is near) a threshold.      
   Should usually be about TSIL_EPSILON, except in very extreme cases. */

#define REDOANALYTIC YES
/* If YES (NO), program does (does not) look for analytic subcases of
   B,S,T,U,V, when the master integral M cannot be done in terms of
   polylogarithms. Should be set to YES, except for testing
   purposes. */

#endif /* TSIL_PARAMS_H */
