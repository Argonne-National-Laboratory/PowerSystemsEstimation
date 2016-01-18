

#ifndef COMMON_INPEST
#define COMMON_INPEST

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

typedef struct {
  DM          dmgen, dmnet; /* DMs to manage generator and network subsystem */
  DM          dmpgrid; /* Composite DM to manage the entire power grid */
  Mat         Ybus; /* Network admittance matrix */
  Vec         V0;  /* Initial voltage vector (Power flow solution) */
  PetscReal   t0,tdisturb,trestore,tfinal; /* inertia initial, disturbance, restore and final times */
  /*PetscScalar Rfault;*/
  PetscReal   dt; /* Initial time and integration time step (or initial time step) */
  PetscInt    neqs_gen,neqs_net,neqs_pgrid;
  PetscBool   saveSol;
  PetscBool   outputCov; /* output prior and data covariances and exit */
  Mat         Sol; /* Matrix to save solution at each time step */
  PetscInt    stepnum;
  //PetscBool   alg_flg;
  PetscReal   t;
  IS          is_diff; /* indices for differential equations */
  IS          is_alg; /* indices for algebraic equations */
  PetscReal   freq_u,freq_l; /* upper and lower frequency limit */
  PetscInt    pow; /* power coefficient used in the cost function */
  Vec         vec_q;
  Vec         proj; /* projects state to observations */
  Mat         obs; /* observations: columns keep time observations, rows keep bus observations */
  char        loadObsFile[1024]; /* path to the file containing observations */
  PetscReal   data_dt; /* time step for data */
  Vec         X0_disturb; /* initial conditions for the system with disturbed inertia */

  PetscReal   misfit; /* accumulate misfit in this term */
  PetscReal   *data_stddev;    /* the noise standard deviation */
  PetscReal   data_noise;    /* data noise at level 'noise' percent */

  PetscReal   prior;    /* accumulate reg in this term */  
  PetscReal   prior_noise; /* the amount of level in prior (controlled by user) */
  PetscReal   *prior_mean;  /* the prior mean */
  PetscReal   *prior_stddev; /* prior standard deviation */

  PetscReal   cost;    /* accumulate cost in this term */
  Mat         J,Jacp; /* Noemi: describe these vars */
  void        *buf1;
} Userctx;


/* Common parameters and --- global, mostly optimization variables --- */

#define freq 60
#define w_s (2*PETSC_PI*freq)

/* Sizes and indices */
extern const PetscInt nbus; /* Number of network buses */
extern const PetscInt ngen; /* Number of generators */
extern const PetscInt nload; /* Number of loads */
extern const PetscInt gbus[3]; /* Buses at which generators are incident */
extern const PetscInt lbus[3]; /* Buses at which loads are incident */

/* Generator real and reactive powers (found via loadflow) */
extern const PetscScalar PG[3];
/* PetscScalar PG[3] = {0.716786142395021,1.630000000000000,0.850000000000000};*/
extern const PetscScalar QG[3];
/* Generator constants */
extern const PetscScalar H0[3];   /* Inertia constant (const declare)*/
extern       PetscScalar H[3];   /* Inertia constant */ 
extern const PetscScalar Rs[3];  /* Stator Resistance */
extern const PetscScalar Xd[3];   /* d-axis reactance */
extern const PetscScalar Xdp[3];  /* d-axis transient reactance */
extern const PetscScalar Xq[3];  /* q-axis reactance Xq(1) set to 0.4360, value given in text 0.0969 */
extern const PetscScalar Xqp[3];  /* q-axis transient reactance */
extern const PetscScalar Td0p[3];  /* d-axis open circuit time constant */
extern const PetscScalar Tq0p[3]; /* q-axis open circuit time constant */
extern PetscScalar M[3]; /* M = 2*H/w_s */
extern PetscScalar D[3]; /* D = 0.1*M */

extern PetscScalar TM[3]; /* Mechanical Torque */
/* Exciter system constants */
extern const PetscScalar KA[3];  /* Voltage regulartor gain constant */
extern const PetscScalar TA[3];     /* Voltage regulator time constant */
extern const PetscScalar KE[3];     /* Exciter gain constant */
extern const PetscScalar TE[3]; /* Exciter time constant */
extern const PetscScalar KF[3];  /* Feedback stabilizer gain constant */
extern const PetscScalar TF[3];    /* Feedback stabilizer time constant */
extern const PetscScalar k1[3];
extern const PetscScalar k2[3];  /* k1 and k2 for calculating the saturation function SE = k1*exp(k2*Efd) */

extern PetscScalar Vref[3];
/* Load constants
  We use a composite load model that describes the load and reactive powers at each time instant as follows
  P(t) = \sum\limits_{i=0}^ld_nsegsp \ld_alphap_i*P_D0(\frac{V_m(t)}{V_m0})^\ld_betap_i
  Q(t) = \sum\limits_{i=0}^ld_nsegsq \ld_alphaq_i*Q_D0(\frac{V_m(t)}{V_m0})^\ld_betaq_i
  where
    ld_nsegsp,ld_nsegsq - Number of individual load models for real and reactive power loads
    ld_alphap,ld_alphap - Percentage contribution (weights) or loads
    P_D0                - Real power load
    Q_D0                - Reactive power load
    V_m(t)              - Voltage magnitude at time t
    V_m0                - Voltage magnitude at t = 0
    ld_betap, ld_betaq  - exponents describing the load model for real and reactive part

    Note: All loads have the same characteristic currently.
*/

extern const PetscScalar PD0_ref[3];         /* const declare corresponding to steady state */
extern       PetscScalar PD0_disturb[3];     /* const declare corresponding to disturbed state */
extern       PetscScalar PD0[3];             /* global variable used in the DAE */
extern const PetscScalar QD0[3];
extern const PetscInt    ld_nsegsp[3];
extern const PetscScalar ld_alphap[3];
extern const PetscScalar ld_betap[3]; 
extern const PetscInt    ld_nsegsq[3];
extern const PetscScalar ld_alphaq[3];
extern const PetscScalar ld_betaq[3]; 


PetscErrorCode ModelSetup(Userctx *user);

/* Integrates the DAE for the given value of P
   Output is the state at final time.
   Additional outputing can be done using PostStepFunction.
 */
PetscErrorCode ForwardSolve(Vec P,void *ctx0, PetscErrorCode (*PostStepFunction)(TS));


PetscErrorCode PreallocateJacobian(Mat J, Userctx *user);
PetscErrorCode SetInitialGuess(Vec X,Userctx *user);
/* \dot{x} - f(x,y)
     g(x,y) = 0
 */
PetscErrorCode IFunction(TS ts,PetscReal t, Vec X, Vec Xdot, Vec F, Userctx *user);

/* Computes F = [-f(x,y);g(x,y)] */
PetscErrorCode ResidualFunction(SNES snes,Vec X, Vec F, Userctx *user);
/* This function is used for solving the algebraic system only .
 It computes the entire F and then zeros out the part corresponding to
   differential equations
 F = [0;g(y)];
*/
PetscErrorCode AlgFunction(SNES snes, Vec X, Vec F, void *ctx);

/*
   J = [a*I-df_dx, -df_dy
        dg_dx, dg_dy]
*/
PetscErrorCode IJacobian(TS ts,PetscReal t,Vec X,Vec Xdot,PetscReal a,Mat A,Mat B,Userctx *user);

/*
   J = [-df_dx, -df_dy
        dg_dx, dg_dy]
*/
PetscErrorCode ResidualJacobian(SNES snes,Vec X,Mat J,Mat B,void *ctx);
PetscErrorCode AlgJacobian(SNES snes,Vec X,Mat A,Mat B,void *ctx);

/* this is a utility to compute a standard normal random variable using the
  * Box-Muller method */
double nrand ();
#endif
