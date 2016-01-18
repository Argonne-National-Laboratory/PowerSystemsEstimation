/*
   Date: 1/17/2016.
   Authors: 
     Noemi Petra (npetra@ucmerced.edu)
     Cosmin G. Petra (petra@mcs.anl.gov
     Emil M. Constantinescu (emconsta@mcs.anl.gov)
     Credit also goes to PETSc 'ex9busopt' examples.

   DISCLAIMER

   THE SOFTWARE IS SUPPLIED .AS IS. WITHOUT WARRANTY OF ANY KIND.

   NEITHER THE UNTED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR UCHICAGO ARGONNE, LLC, NOR ANY OF THEIR EMPLOYEES, NOR ANY OF THE AUTOHRS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

static char help[] = "Power grid inertia estimation of WECC 9 bus system.\n\
This example is based on the 9-bus (node) example given in the book Power\n\
Systems Dynamics and Stability (Chapter 7) by P. Sauer and M. A. Pai. and \n\
'ex9busopt' example from PETSc distribution.                              \n\
The power grid in this example consists of 9 buses (nodes), 3 generators,\n\
3 loads, and 9 transmission lines. The network equations are written\n\
in current balance form using rectangular coordinates.\n\n";


/*
   The equations for the stability analysis are described by the DAE

   \dot{x} = f(x,y,t)
     0     = g(x,y,t)

   where the generators are described by differential equations, while the algebraic
   constraints define the network equations.

   The generators are modeled with a 4th order differential equation describing the electrical
   and mechanical dynamics. Each generator also has an exciter system modeled by 3rd order
   diff. eqns. describing the exciter, voltage regulator, and the feedback stabilizer
   mechanism.

   The network equations are described by nodal current balance equations.
    I(x,y) - Y*V = 0

   where:
    I(x,y) is the current injected from generators and loads.
      Y    is the admittance matrix, and
      V    is the voltage vector
*/

/*
   Include "petscts.h" so that we can use TS solvers.  Note that this
   file automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscksp.h   - linear solvers
*/


#include "common-inpest.h"
#include <petsctao.h>
#include <time.h>
//#include <petsctime.h>

PetscErrorCode FormFunction(Tao,Vec,PetscReal*,void*);
PetscErrorCode FormFunctionGradient(Tao,Vec,PetscReal*,Vec,void*);
PetscErrorCode InitializeData(const PetscScalar* P, void *ctx0, double noise, PetscScalar data_dt);
//PetscErrorCode PetscGetTime(PetscLogDouble *t);

/* ********************************************************************** */
/* Noemi: describe these functions in detail! */
/* ********************************************************************** */

/* Matrix JacobianP is constant so that it only needs to be evaluated once */
#undef __FUNCT__
#define __FUNCT__ "RHSJacobianP"
static PetscErrorCode RHSJacobianP(TS ts,PetscReal t,Vec X,Mat A,void *ctx0)
{
  PetscErrorCode ierr;
  PetscScalar    a;
  PetscInt       row,col;
  Userctx        *ctx=(Userctx*)ctx0;
  PetscInt       idx=0;
  Vec            Xgen,Xnet;
  PetscScalar    *xgen,*xnet;

  PetscScalar    Eqp,Edp;
  PetscScalar    Id,Iq;

  PetscFunctionBeginUser;

  /*if (ctx->jacp_flg) { delete me */
  ierr = MatZeroEntries(A);CHKERRQ(ierr);

    //recompute since this changes
  M[0] = 2*H[0]/w_s; M[1] = 2*H[1]/w_s; M[2] = 2*H[2]/w_s;
  D[0] = 0.1*M[0]; D[1] = 0.1*M[1]; D[2] = 0.1*M[2];

  ierr = DMCompositeGetLocalVectors(ctx->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  ierr = DMCompositeScatter(ctx->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);

  /* Generator subsystem initialization */
  ierr = VecGetArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);

  for (col=0;col<ngen;col++) {

    Eqp   = xgen[idx];
    Edp   = xgen[idx+1];
    Id    = xgen[idx+4];
    Iq    = xgen[idx+5];
    TM[col] = PG[col];

    a = -1.0 *(TM[col] - Edp*Id - Eqp*Iq - (Xqp[col]-Xdp[col])*Id*Iq)/M[col]/H[col];

    row  = 9*col+3; // E is in the 4th equation hence +3
    idx = idx + 9;  //9 equations per generator, hence move to next gen

    ierr = MatSetValues(A,1,&row,1,&col,&a,INSERT_VALUES);CHKERRQ(ierr);
 }

  /* ctx->jacp_flg = PETSC_FALSE; delete me */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xgen,&xgen);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);

  ierr = DMCompositeGather(ctx->dmpgrid,X,INSERT_VALUES,Xgen,Xnet);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(ctx->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  //MatView(A,PETSC_VIEWER_STDOUT_SELF);

  /* }delete me */
  PetscFunctionReturn(0);
}

/* This function implements the cost  */
/*#undef __FUNCT__
#define __FUNCT__ "CostIntegrand"
static PetscErrorCode CostIntegrand(TS ts,PetscReal t,Vec U,Vec R,Userctx *user)
{
  PetscErrorCode    ierr;
  PetscScalar       *r;
  const PetscScalar *u;
  Vec               Xgen,Xnet;
  PetscScalar       *xnet,*proj_vec;
  PetscInt          obs_len,i,idx;
  PetscReal         step_num;
  PetscScalar       *mat;

  PetscFunctionBegin;
  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->dmpgrid,U,Xgen,Xnet);CHKERRQ(ierr);

  ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(R,&r);CHKERRQ(ierr);
  r[0] = 0.;

  t = t - user->tdisturb;
  step_num = round(t / user->data_dt);
  if(fabs(step_num - t/user->data_dt)<= 1e-6*user->dt) {
    ierr     = VecGetArray(user->proj,&proj_vec);CHKERRQ(ierr);
    ierr     = VecGetSize(user->proj, &obs_len);CHKERRQ(ierr);
    idx      = 2*obs_len*(PetscInt) step_num;

    //ierr     = TSGetSolution(ts,&U);CHKERRQ(ierr);

    ierr     = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
    ierr     = DMCompositeScatter(user->dmpgrid,U,Xgen,Xnet);CHKERRQ(ierr);
    ierr     = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
    ierr     = MatDenseGetArray(user->obs,&mat);CHKERRQ(ierr);

    for(i=0;i<obs_len; i++) {
      r[0] += 0.5*pow(xnet[2*((int)proj_vec[i])]   - mat[idx+2*i],  2);  // obs[2*i, num_obs]
      r[0] += 0.5*pow(xnet[2*((int)proj_vec[i])+1] - mat[idx+2*i+1],2);  // obs[2*i+1, num_obs]
    }
  }
  //printf("CostIntegr(%6.3f)=%g\n", t, r[0]);

  ierr     = MatDenseRestoreArray(user->obs,&mat);CHKERRQ(ierr);
  ierr     = VecRestoreArray(user->proj,&proj_vec);CHKERRQ(ierr);

  ierr     = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);
  ierr     = VecRestoreArray(R,&r);CHKERRQ(ierr);
  ierr     = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr     = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
*/

/* This function implements the  d (mifit) / dy, where y = (x,y) */
#undef __FUNCT__
#define __FUNCT__ "DRDYFunction"
static PetscErrorCode DRDYFunction(TS ts,PetscReal t,Vec U,Vec *drdy,Userctx *user)
{
  PetscErrorCode ierr;
  /*
  Vec            Xgen,Xnet,Dgen,Dnet;
  PetscScalar    *xnet,*dgen,*proj_vec;
  PetscInt       obs_len,i,idx;
  PetscReal      step_num;
  PetscScalar    *mat;
  */
  PetscFunctionBegin;

  ierr = VecZeroEntries(drdy[0]);CHKERRQ(ierr);
  /*  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Dgen,&Dnet);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->dmpgrid,U,Xgen,Xnet);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->dmpgrid,drdy[0],Dgen,Dnet);CHKERRQ(ierr);

  ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
  ierr = VecGetArray(Dgen,&dgen);CHKERRQ(ierr);

  t = t - user->tdisturb;
  step_num = round(t / user->data_dt);
  if(fabs(step_num - t/user->data_dt)<= 1e-6*user->dt) {
    ierr     = VecGetArray(user->proj,&proj_vec);CHKERRQ(ierr);
    ierr     = VecGetSize(user->proj, &obs_len);CHKERRQ(ierr);
    idx      = 2*obs_len*(PetscInt) step_num;

    ierr     = TSGetSolution(ts,&U);CHKERRQ(ierr);

    ierr     = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
    ierr     = DMCompositeScatter(user->dmpgrid,U,Xgen,Xnet);CHKERRQ(ierr);
    ierr     = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
    ierr     = MatDenseGetArray(user->obs,&mat);CHKERRQ(ierr);

    for(i=0;i<obs_len; i++) {
      //printf("B^T (B x -d): t=%g index=%d proj_len=%d proj[i] = %d \n", t, idx, obs_len, (int)proj_vec[i]);
      dgen[2*((int)proj_vec[i])] += xnet[2*((int)proj_vec[i])]   - mat[idx+2*i];    // obs[2*i, num_obs]
      dgen[2*((int)proj_vec[i])] += xnet[2*((int)proj_vec[i])+1] - mat[idx+2*i+1];  // obs[2*i+1, num_obs]
    }
  }

  ierr     = MatDenseRestoreArray(user->obs,&mat);CHKERRQ(ierr);
  ierr     = VecRestoreArray(user->proj,&proj_vec);CHKERRQ(ierr);

  ierr = VecRestoreArray(Dgen,&dgen);CHKERRQ(ierr);
  ierr = DMCompositeGather(user->dmpgrid,drdy[0],INSERT_VALUES,Dgen,Dnet);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Dgen,&Dnet);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  */
  PetscFunctionReturn(0);
}

/* This function adds d (misfit)(t)/dy to drdy, where y = (x,y) */
#undef __FUNCT__
#define __FUNCT__ "AddDRDY"
static PetscErrorCode AddDRDY(PetscReal t,Vec X,Vec *drdy,Userctx *user)
{
  PetscErrorCode ierr;
  Vec            Xgen,Xnet,Dgen,Dnet;
  PetscScalar    *xnet,*dnet,*proj_vec;
  PetscInt       obs_len,i,idx;
  PetscReal      step_num, aux1, aux2;
  PetscScalar    *mat;

  PetscFunctionBegin;

  //printf("AddDRDY lambda\n");
  //ierr = VecView(*drdy, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Dgen,&Dnet);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->dmpgrid,*drdy,Dgen,Dnet);CHKERRQ(ierr);

  ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
  ierr = VecGetArray(Dnet,&dnet);CHKERRQ(ierr);

  ierr     = VecGetArray(user->proj,&proj_vec);CHKERRQ(ierr);
  ierr     = VecGetSize(user->proj, &obs_len);CHKERRQ(ierr);

  //!step_num = round((t-user->tdisturb) / user->data_dt);
  step_num = round((t-user->trestore) / user->data_dt);
  idx      = 2*obs_len* (PetscInt)step_num;  
  ierr     = MatDenseGetArray(user->obs,&mat);CHKERRQ(ierr);  


  for(i=0;i<obs_len; i++) {
    //printf("B^Tnoisecovinv(B x -d):t=%g index=%d proj_len=%d proj[i] = %d \n", t, idx, obs_len, (int)proj_vec[i]);
    aux1 = (xnet[2*((int)proj_vec[i])]   - mat[idx+2*i])
      / pow(user->data_stddev[i],2);    // obs[2*i, num_obs]
    dnet[2*((int)proj_vec[i])] += aux1;
    
    aux2 = (xnet[2*((int)proj_vec[i])+1] - mat[idx+2*i+1])
      / pow(user->data_stddev[i+1],2);// obs[2*i+1, num_obs]
    dnet[2*((int)proj_vec[i])+1] += aux2;

    //printf("IC Adjoint: t=%g [%g] [%g]\n", t, aux1, aux2);
  }

  ierr     = MatDenseRestoreArray(user->obs,&mat);CHKERRQ(ierr);
  ierr     = VecRestoreArray(user->proj,&proj_vec);CHKERRQ(ierr);

  ierr = VecRestoreArray(Dnet,&dnet);CHKERRQ(ierr);
  ierr = DMCompositeGather(user->dmpgrid,*drdy,INSERT_VALUES,Dgen,Dnet);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Dgen,&Dnet);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);

  //printf("AddDRDY lambda at the end\n");
  //ierr = VecView(*drdy, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* This function implements the  d (misfit) / dp */
#undef __FUNCT__
#define __FUNCT__ "DRDPFunction"
static PetscErrorCode DRDPFunction(TS ts,PetscReal t,Vec U,Vec *drdp,Userctx *user)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = VecZeroEntries(drdp[0]);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/* ******************************************************************************* */

/*#undef __FUNCT__
#define __FUNCT__ "MonitorUpdateQ"
static PetscErrorCode MonitorUpdateQ(TS ts,PetscInt stepnum,PetscReal time,Vec X,void *ctx0)
{
  PetscErrorCode ierr;
  Vec            C,*Y;
  PetscInt       Nr;
  PetscReal      h,theta; 
  Userctx        *ctx=(Userctx*)ctx0;
 
  PetscFunctionBegin;
  theta = 0.5;
  ierr = TSGetStages(ts,&Nr,&Y);CHKERRQ(ierr);
  ierr = TSGetTimeStep(ts,&h);CHKERRQ(ierr);
  ierr = VecDuplicate(ctx->vec_q,&C);CHKERRQ(ierr);
  / * compute integrals * /
  if (stepnum>0) {
    ierr = CostIntegrand(ts,time,X,C,ctx);CHKERRQ(ierr);
    ierr = VecAXPY(ctx->vec_q,h*theta,C);CHKERRQ(ierr);
    ierr = CostIntegrand(ts,time+h*theta,Y[0],C,ctx);CHKERRQ(ierr);
    ierr = VecAXPY(ctx->vec_q,h*(1-theta),C);CHKERRQ(ierr);  
  }
  ierr = VecDestroy(&C);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
*/

/* Saves the solution at each time to a matrix */
#undef __FUNCT__
#define __FUNCT__ "SetSolution"
PetscErrorCode SetSolution(Userctx* user, PetscReal t, Vec X)
{
  PetscErrorCode ierr;
  PetscScalar    *xgen,*xnet,*mat;
  PetscInt       idx, len_xgen, len_xnet;
  Vec Xgen, Xnet;

  PetscFunctionBegin;

  if(user->saveSol) {

    //printf("SaveSol at t=%g\n", t);

    ierr = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
    ierr = DMCompositeScatter(user->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);
    ierr = VecGetSize(Xgen, &len_xgen);CHKERRQ(ierr);
    ierr = VecGetSize(Xnet, &len_xnet);CHKERRQ(ierr);
    ierr = VecGetArray(Xgen,&xgen);CHKERRQ(ierr);
    ierr = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);

    idx      = ((PetscInt) round(t / user->dt))*(user->neqs_pgrid+1);
    ierr     = MatDenseGetArray(user->Sol,&mat);CHKERRQ(ierr);
    mat[idx] = t;


    ierr     = PetscMemcpy(mat+idx+1,         xgen,len_xgen*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr     = PetscMemcpy(mat+idx+1+len_xgen,xnet,len_xnet*sizeof(PetscScalar));CHKERRQ(ierr);


    ierr     = MatDenseRestoreArray(user->Sol,&mat);CHKERRQ(ierr);


    ierr = VecRestoreArray(Xgen,&xgen);CHKERRQ(ierr);
    ierr = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);
    ierr = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SetObservation"
PetscErrorCode SetObservation(Userctx* user, PetscReal t, Vec X)
{
  PetscErrorCode    ierr;
  PetscInt          idx,obs_len,i;
  PetscScalar       *mat;
  Vec               Xgen,Xnet;
  PetscScalar       *xnet;//, *xgen;
  const PetscScalar *x;
  PetscScalar       *proj_vec;
  PetscReal         step_num;

  PetscFunctionBegin;

  /* observations are generated only after restore */
  /* t = t - user->tdisturb;*/
  t = t - user->trestore;
  if(t >= -1e-6*user->dt) {
    step_num = round(t / user->data_dt);
    if(fabs(step_num - t/user->data_dt)<= 1e-6*user->dt) {
      
      idx      = (PetscInt) step_num;
      ierr     = VecGetSize(user->proj, &obs_len);CHKERRQ(ierr);
      
      //printf("SetObservation: t=%g index=%d proj_len=%d \n", t, idx, obs_len); 
      
      idx      = idx*2*obs_len;
      
      ierr     = MatDenseGetArray(user->obs,&mat);CHKERRQ(ierr);
      ierr     = VecGetArrayRead(X,&x);CHKERRQ(ierr);
      ierr     = VecGetArray(user->proj, &proj_vec);CHKERRQ(ierr);
      
      ierr     = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
      ierr     = DMCompositeScatter(user->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);
      ierr     = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
      
      //get pointer to xnet and to user->proj
      for(i=0;i<obs_len; i++) {
	//printf("proj[%d]=%d\n", i, (int)proj_vec[i]);
	*(mat+idx+2*i)   = xnet[2*((int)proj_vec[i])];
	*(mat+idx+2*i+1) = xnet[2*((int)proj_vec[i])+1];
      }
      //printf("%g \n", xnet[0]);
      ierr     = MatDenseRestoreArray(user->obs,&mat);CHKERRQ(ierr);
      ierr     = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
      ierr     = VecRestoreArray(user->proj,&proj_vec);CHKERRQ(ierr);
      ierr     = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);
      ierr     = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}


/* Saves the solution at each time to a matrix */
#undef __FUNCT__
#define __FUNCT__ "SaveObservation"
PetscErrorCode SaveObservation(TS ts)
{
  PetscErrorCode    ierr;
  Userctx           *user;
  Vec               X;
  PetscReal         t;


  PetscFunctionBegin;

  ierr     = TSGetApplicationContext(ts,&user);CHKERRQ(ierr);
  ierr     = TSGetTime(ts,&t);CHKERRQ(ierr);
  ierr     = TSGetSolution(ts,&X);CHKERRQ(ierr);

  ierr     = SetObservation(user,t,X);CHKERRQ(ierr);
  ierr     = SetSolution(user,t,X);CHKERRQ(ierr);

  user->stepnum++;
  PetscFunctionReturn(0);
}


/* Evaluate the misfit at each observation time */
#undef __FUNCT__
#define __FUNCT__ "EvalMisfit"
PetscErrorCode EvalMisfit(TS ts)
{
  PetscErrorCode    ierr;
  Userctx           *user;
  Vec               X;
  PetscReal         t;
  Vec               Xgen,Xnet;
  PetscScalar       *xnet, *proj_vec;
  PetscInt          obs_len,i,idx;
  PetscReal         step_num;
  PetscScalar       *mat;

  PetscFunctionBegin;

  ierr     = TSGetApplicationContext(ts,&user);CHKERRQ(ierr);
  ierr     = TSGetTime(ts,&t);CHKERRQ(ierr);

  ierr     = TSGetSolution(ts,&X);CHKERRQ(ierr);


  //!t = t - user->tdisturb;
  t = t - user->trestore;
  step_num = round(t / user->data_dt);
  if(fabs(step_num - t/user->data_dt)<= 1e-6*user->dt) {
 
    ierr     = VecGetArray(user->proj,&proj_vec);CHKERRQ(ierr);
    ierr     = VecGetSize(user->proj, &obs_len);CHKERRQ(ierr);
    idx      = 2*obs_len*(PetscInt) step_num;

    ierr     = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
    ierr     = DMCompositeScatter(user->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);
    ierr     = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
    ierr     = MatDenseGetArray(user->obs,&mat);CHKERRQ(ierr);

    for(i=0;i<obs_len; i++) {
      user->misfit += 0.5*pow(xnet[2*((int)proj_vec[i])]   - mat[idx+2*i],  2) / pow(user->data_stddev[i], 2); // obs[2*i, num_obs]
      user->misfit += 0.5*pow(xnet[2*((int)proj_vec[i])+1] - mat[idx+2*i+1],2)
	/ pow(user->data_stddev[i+1], 2);// obs[2*i+1, num_obs]
    }

    //printf("evm %g %g idx=%d\n", xnet[0], mat[idx], idx);
    ierr     = MatDenseRestoreArray(user->obs,&mat);CHKERRQ(ierr);
    ierr     = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);
    ierr     = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
    ierr     = VecRestoreArray(user->proj,&proj_vec);CHKERRQ(ierr);
  }

  user->stepnum++;
  PetscFunctionReturn(0);
}


/* Evaluate the regularization term at current parameter */
#undef __FUNCT__
#define __FUNCT__ "EvalReg"
PetscErrorCode EvalReg(Userctx *user, Vec P)
{
  PetscErrorCode    ierr;
  PetscScalar       *x_ptr;

  PetscFunctionBegin;

  ierr  = VecGetArray(P,&x_ptr);CHKERRQ(ierr);

  user->prior = 0.5*(pow(x_ptr[0]-user->prior_mean[0],2) / pow(user->prior_stddev[0],2)+
		   pow(x_ptr[1]-user->prior_mean[1],2) / pow(user->prior_stddev[1],2)+
		   pow(x_ptr[2]-user->prior_mean[2],2) / pow(user->prior_stddev[2],2));
  ierr  = VecRestoreArray(P,&x_ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Evaluate the regularization term at current parameter */
#undef __FUNCT__
#define __FUNCT__ "AddRegGradient"
PetscErrorCode AddRegGradient(Userctx *user, Vec P, Vec grad)
{
  PetscErrorCode    ierr;
  PetscScalar       *p, *g;

  PetscFunctionBegin;

  ierr  = VecGetArray(P,&p);CHKERRQ(ierr);
  ierr  = VecGetArray(grad,&g);CHKERRQ(ierr);

  g[0] += (p[0]-user->prior_mean[0]) / pow(user->prior_stddev[0],2);
  g[1] += (p[1]-user->prior_mean[1]) / pow(user->prior_stddev[1],2);
  g[2] += (p[2]-user->prior_mean[2]) / pow(user->prior_stddev[2],2);

  ierr  = VecRestoreArray(P,&p);CHKERRQ(ierr);
  ierr  = VecRestoreArray(grad,&g);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "ComputeSensiP"
PetscErrorCode ComputeSensiP(Vec lambda,Vec mu,Vec *DICDP,Userctx *user)
{
  PetscErrorCode ierr;
  PetscScalar    *x,*y,sensip;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = VecGetArray(lambda,&x);CHKERRQ(ierr);
  ierr = VecGetArray(mu,&y);CHKERRQ(ierr);

  for (i=0;i<3;i++) {
    ierr   = VecDot(lambda,DICDP[i],&sensip);CHKERRQ(ierr);
    sensip = sensip+y[i];
    ierr   = PetscPrintf(PETSC_COMM_WORLD,"\n sensitivity wrt %D th parameter: %g \n",
			 i,(double)sensip);CHKERRQ(ierr);
     y[i] = sensip;
  }
  ierr = VecRestoreArray(mu,&y);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateHessianFD"
PetscErrorCode EvaluateHessianFD(Tao tao, Vec p, Userctx* user)
{
  PetscErrorCode ierr;
  PetscReal f, expo, *eps_arr;
  Vec G, G_eps;
  int i,j;
  PetscReal minusOne=-1., one=1.;
  Vec p_eps;
  PetscScalar test_eps;

  FILE* fOut = fopen("Hessian.txt", "a+");
  if(fOut==NULL)
    PetscFunctionReturn(-1);

  //alocate the vector for gradients: at current point and at perturbed point
  ierr = VecDuplicate(p, &G);CHKERRQ(ierr);
  ierr = VecDuplicate(p, &G_eps);CHKERRQ(ierr);

  //gradient at current point
  ierr = FormFunctionGradient(tao, p, &f, G, user);

  //alocate the perturbed x
  ierr = VecDuplicate(p, &p_eps);CHKERRQ(ierr);

  /* testing flag: 0 means fix epsilon; 1 means try several epsilons*/
  test_eps = 0;
  for(i=0;i<3;i++) {
    printf("column %d\n", i);
    //test different expo to see what eps works better
    for(expo=1e-2; expo>1e-8; expo/=3) {
      if (test_eps == 0){
	expo = 1.0e-4;
      }
      ierr =  VecCopy(p, p_eps);
      ierr = VecGetArray(p_eps, &eps_arr);CHKERRQ(ierr);
      eps_arr[i] += expo;
      ierr = VecRestoreArray(p_eps, &eps_arr);CHKERRQ(ierr);

      //gradient at perturbed point
      ierr = FormFunctionGradient(tao, p_eps, &f, G_eps, user);

      //i-th column of the Hessian
      //G_eps = G_eps-G
      ierr = VecAXPBY(G_eps, minusOne, one, G);CHKERRQ(ierr);
      //G_eps = G_eps/eps
      ierr = VecScale(G_eps, 1/expo);CHKERRQ(ierr);

      printf("eps=%7.4e ", expo);
      ierr = VecGetArray(G_eps,  &eps_arr);CHKERRQ(ierr);
      for(j=0; j<3; j++){
	if (test_eps == 1){
	  printf( "%18.12f   ", eps_arr[j]);
	}
	else{
	  if (j==0){
	      fprintf(fOut, " %.12f %.12f %.12f %18.12f ", user->tfinal, user->data_noise, PD0_disturb[0], eps_arr[j]);
	      printf( " %.12f %.12f %.12f %18.12f ", user->tfinal, user->data_noise, PD0_disturb[0], eps_arr[j]);
	    }
	  else{
	    fprintf(fOut, " %18.12f ", eps_arr[j]);
	    printf( " %18.12f ", eps_arr[j]);
	  }
	}
      }
      if (test_eps == 1){
	printf("\n");
      }
      else{
	fprintf(fOut, "\n");
	printf("\n");
	break;
      }
      ierr = VecRestoreArray(G_eps, &eps_arr);CHKERRQ(ierr);
    }
  }
  fclose(fOut);
  //destroy vectors

  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&G_eps);CHKERRQ(ierr);
  ierr = VecDestroy(&p_eps);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Userctx            user;
  Vec                p;
  PetscScalar        *x_ptr;
  PetscErrorCode     ierr;
  PetscMPIInt        size;
  PetscInt           i,numDataBuses;
  KSP                ksp;
  PC                 pc;
  Tao                tao;
  TaoConvergedReason reason;
  Vec                lowerb,upperb;
  PetscViewer        viewer;
  PetscScalar   *proj_vec;

  //PetscLogDouble     t0,t1;

  /* time the inversion process */
  //ierr = PetscGetTime(&t0);CHKERRQ(ierr);

  ierr = PetscInitialize(&argc,&argv,"petscoptions",help);CHKERRQ(ierr);
  PetscFunctionBeginUser;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only for sequential runs");

  ierr = ModelSetup(&user);CHKERRQ(ierr);

  /* hard code the data projection here - for now assume data at all buses */
  ierr = VecCreateSeq(PETSC_COMM_WORLD,nbus,&user.proj);CHKERRQ(ierr);
  /*ierr = VecCreateSeq(PETSC_COMM_WORLD,4,&user.proj);CHKERRQ(ierr);*/
  ierr = VecGetArray(user.proj,&proj_vec);CHKERRQ(ierr);
  for(i=0; i<nbus; i++) {
    proj_vec[i]=i;
  }

  srand( time(NULL) + rand () ); 

  //VecView(user.proj, PETSC_VIEWER_STDOUT_WORLD);
  /* -- 2 5 6 8 */
  /* -- proj_vec[0]=1; proj_vec[1]=4; proj_vec[2]=5; proj_vec[3]=7; */
  ierr = VecRestoreArray(user.proj,&proj_vec);CHKERRQ(ierr);

  /* allocate/set the prior mean and its standard deviation */
  ierr = PetscMalloc(3*sizeof(PetscScalar), &user.prior_mean);
  ierr = PetscMalloc(3*sizeof(PetscScalar), &user.prior_stddev);


  /*{23.64,6.4,3.01};*/
  user.prior_mean[0] = 24.0;
  user.prior_mean[1] = 6.0;
  user.prior_mean[2] = 3.1;
  for(i=0; i<3; i++) user.prior_stddev[i] = user.prior_mean[i]*user.prior_noise;

  /* Create matrix to store solution */
  if(user.saveSol) {
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,
			     user.neqs_pgrid+1, 
			     (PetscInt) round((user.tfinal-user.t0)/user.dt+1),
			     NULL,
			     &user.Sol);
    CHKERRQ(ierr);
  }
  printf("Num cols=%d\n", (PetscInt) round((user.tfinal-user.t0)/user.dt+1));

  /* *********************************
   * Generate/load observations
   **********************************/
  ierr = VecGetSize(user.proj, &numDataBuses);CHKERRQ(ierr);
  /* Create matrix to save solutions at each time step */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, 
			   2*numDataBuses, 
			   //(PetscInt) round((user.tfinal-user.tdisturb)/user.data_dt)+1, 
			   (PetscInt) round((user.tfinal-user.trestore)/user.data_dt)+1, 
			   NULL, 
			   &user.obs);
  CHKERRQ(ierr);

  ierr = InitializeData(H0, &user, user.data_noise, user.data_dt);CHKERRQ(ierr);

  if(0==strlen(user.loadObsFile)) {
    /*  save observations */
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"obs-perturbed.bin",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = MatView(user.obs,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    printf("Observations generated.\n");
  }

  if(user.saveSol) {
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"out_pert.bin",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = MatView(user.Sol,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = MatDestroy(&user.Sol);CHKERRQ(ierr);

    CHKERRQ(ierr);
  }

  if(user.outputCov) {
    printf("The diagonal of the data noise covariance matrix (%g absolute noise) is:\n", user.data_noise);
    for(i=0; i<2*numDataBuses; i++) printf("%18.12f ", user.data_stddev[i]*user.data_stddev[i]); printf("\n");

    printf("The prior mean is: ");
    for(i=0; i<3; i++) printf("%18.12f ", user.prior_mean[i]); printf("\n");

    printf("The diagonal of the prior covariance matrix (%g relative noise) is:\n", user.prior_noise);
    for(i=0; i<3; i++) printf("%18.12f ", user.prior_stddev[i]*user.prior_stddev[i]); printf("\n");
    goto finalize;
  }

  /* ***************************************
   *    Optimization phase
   * ***************************************/
  /* Create TAO solver and set desired solution method */
  ierr = TaoCreate(PETSC_COMM_WORLD,&tao);CHKERRQ(ierr);
  ierr = TaoSetType(tao,TAOBLMVM);CHKERRQ(ierr);
  /*
     Optimization starts
  */

  printf("Starting optimization...\n");

  /* PetscScalar H_disturb[3]= {25.,6.4,3.01};  New inertia (after tdisturb) to be estimated */

  /* Set initial solution guess */
  ierr = VecCreateSeq(PETSC_COMM_WORLD,3,&p);CHKERRQ(ierr);
  ierr = VecGetArray(p,&x_ptr);CHKERRQ(ierr);
  //x_ptr[0] = H0[0]; x_ptr[1] = H0[1]; x_ptr[2] = H0[2];
  x_ptr[0] = H0[0]*1.1; x_ptr[1] = H0[1]*1.1; x_ptr[2] = H0[2]*1.1;
  ierr = VecRestoreArray(p,&x_ptr);CHKERRQ(ierr);

  ierr = TaoSetInitialVector(tao,p);CHKERRQ(ierr);

  /* Set routine for function and gradient evaluation */
  //ierr = TaoSetObjectiveRoutine(tao,FormFunction,(void *)&user);CHKERRQ(ierr);
  //ierr = TaoSetGradientRoutine(tao,TaoDefaultComputeGradient,(void *)&user);CHKERRQ(ierr);

  /* Sets the cost and gradient evaluation routine for minimization */
  ierr = TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,&user);CHKERRQ(ierr);

  /* Set bounds for the optimization */
  ierr = VecDuplicate(p,&lowerb);CHKERRQ(ierr);
  ierr = VecDuplicate(p,&upperb);CHKERRQ(ierr);
  ierr = VecGetArray(lowerb,&x_ptr);CHKERRQ(ierr);
  x_ptr[0] = 20.64; x_ptr[1] = 5.4; x_ptr[2] = 2.01;
  ierr = VecRestoreArray(lowerb,&x_ptr);CHKERRQ(ierr);
  ierr = VecGetArray(upperb,&x_ptr);CHKERRQ(ierr);
  x_ptr[0] = 25.64; x_ptr[1] = 7.4; x_ptr[2] = 4.01;
  ierr = VecRestoreArray(upperb,&x_ptr);CHKERRQ(ierr);
  ierr = TaoSetVariableBounds(tao,lowerb,upperb);

  /* Check for any TAO command line options */
  ierr = TaoSetFromOptions(tao);CHKERRQ(ierr);
  ierr = TaoGetKSP(tao,&ksp);CHKERRQ(ierr);
  if (ksp) {
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
  }

  //ierr = TaoSetTolerances(tao,1e-8,1e-6,1e-8,1e-6,1e-4);
  ierr = TaoSetTolerances(tao,1e-8,1e-8,1e-8,1e-8,1e-6);
  //ierr = TaoSetGradientTolerances(tao,1e-8, 1e-6, 1e-6);

  /* SOLVE the estimation problem */
  ierr = TaoSolve(tao); CHKERRQ(ierr);
  /* Get information on termination */

  printf("--- optimization done\n");
  /* time the inversion process */
  //ierr = PetscGetTime(&t1);CHKERRQ(ierr);
 
  //printf("elapsed_time  %f seconds\n", t1 - t0);

  ierr = TaoGetConvergedReason(tao,&reason);CHKERRQ(ierr);
  if (reason <= 0){
    ierr=PetscPrintf(MPI_COMM_WORLD, "Try another method! \n");CHKERRQ(ierr);
  }

  /*ierr = VecView(p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);*/
  ierr = VecGetArray(p,&x_ptr);CHKERRQ(ierr);
  printf("inertia-out: %.12f %.12f %.12f\n", x_ptr[0], x_ptr[1], x_ptr[2]);
  ierr = VecRestoreArray(p,&x_ptr);CHKERRQ(ierr);

  //ierr =  EvaluateHessianFD(tao, p, &user);CHKERRQ(ierr); 

    /* Free TAO data structures */
  ierr = TaoDestroy(&tao);CHKERRQ(ierr);
  ierr = VecDestroy(&lowerb);CHKERRQ(ierr);
  ierr = VecDestroy(&upperb);CHKERRQ(ierr);

 finalize:
  ierr = MatDestroy(&user.obs);CHKERRQ(ierr);
  ierr = VecDestroy(&user.X0_disturb);CHKERRQ(ierr);
  ierr = PetscFree(user.data_stddev);CHKERRQ(ierr);
  PetscFree(user.prior_mean);
  PetscFree(user.prior_stddev);

  ierr = DMDestroy(&user.dmgen);CHKERRQ(ierr);
  ierr = DMDestroy(&user.dmnet);CHKERRQ(ierr);
  ierr = DMDestroy(&user.dmpgrid);CHKERRQ(ierr);
  ierr = ISDestroy(&user.is_diff);CHKERRQ(ierr);
  ierr = ISDestroy(&user.is_alg);CHKERRQ(ierr);

  ierr = MatDestroy(&user.J);CHKERRQ(ierr);
  ierr = MatDestroy(&user.Jacp);CHKERRQ(ierr);
  ierr = MatDestroy(&user.Ybus);CHKERRQ(ierr);
  ierr = VecDestroy(&user.V0);CHKERRQ(ierr);

  ierr = VecDestroy(&p);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

/* ------------------------------------------------------------------ */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/*
   FormFunction - Evaluates the function and corresponding gradient.

   Input Parameters:
   tao - the Tao context
   X   - the input vector
   ptr - optional user-defined context, as set by TaoSetObjectiveAndGradientRoutine()

   Output Parameters:
   f   - the newly evaluated function
*/
PetscErrorCode FormFunction(Tao tao,Vec P,PetscReal *f,void *ctx0)
{
  PetscErrorCode ierr;
  Userctx        *ctx = (Userctx*)ctx0;
  

  ctx->misfit=0.0;
  ctx->stepnum = 0;

  
  ierr = ForwardSolve(P,ctx0,EvalMisfit);CHKERRQ(ierr);


  /* Compute the cost = misfit + reg*/
  EvalReg(ctx, P);
  *f  = ctx->prior;

  *f += ctx->misfit;
  //printf("misfit=%.6f reg=%.6f\n", ctx->misfit, ctx->reg);

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "InitializeData"
/*
 InitializeData - Generate/Load observations by solving the forward problem with parameters
 P and adding noise to the DAE solution (xnet[i]=xnet[i]*(1+noise*randn))

 Input Parameters:
 P    - the input vector of parameters (inertia of the generators)
 ctx0 - optional user-defined context
 proj - contains the bus indexes at which the data is generated (that is buses where data is available),
 a subset of {0,1,...,nbus}
 noise- ammount of noise level

 Output Parameters:
 obs  - matrix of "observations" (perturbed DAE solution) at 0, data_dt, 2data_dt, ..., ctx->tfinal,
 ith column containing the observations specified by 'proj' at time i (i=0,...,ctx->tfinal/data_dt)
*/

PetscErrorCode InitializeData(const PetscScalar* P, void *ctx0, double noise, PetscScalar data_dt)
{
  TS             ts;
  SNES           snes_alg;
  PetscErrorCode ierr;
  Userctx        *ctx = (Userctx*)ctx0;
  Vec            X;
  Mat            J;
  Vec            F_alg;
  Vec            Xdot;
  PetscInt       i,j,m,n;
  PetscReal      *mat;
  //PetscReal      temp;
  PetscViewer    obsView;

  PetscFunctionBegin;

  H[0] = P[0];
  H[1] = P[1];
  H[2] = P[2];
  printf("InitializeData: x=[%.14f, %.14f, %.14f], obs_noise=%5.3f Nt = %4.2f, Nobs=%4.2f\n",  H[0],  H[1], H[2], noise, ((ctx->tfinal-ctx->t0)/ctx->dt)+1, ((ctx->tfinal-ctx->trestore)/ctx->data_dt)+1);

  if(ctx->t0 > ctx->tdisturb) {
    printf("t0 cannot be greater than tdisturb\n");
    PetscFunctionReturn(-1);
  }
  if( (ctx->tdisturb >= ctx->trestore-1.0e-8) || (ctx->tdisturb >= ctx->tfinal-1.0e-8) ) {
    printf("tdisturb should be less than trestore and tfinal\n");
    PetscFunctionReturn(-1);
  }

  //use the reference PD0 from t0 to t_disturb to ensure steady state
  for(i=0; i<3; i++) PD0[i] = PD0_ref[i];

  ctx->stepnum = 0;

  ierr = DMCreateGlobalVector(ctx->dmpgrid,&X);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
  ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,ctx->neqs_pgrid,ctx->neqs_pgrid);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J);CHKERRQ(ierr);
  ierr = PreallocateJacobian(J,ctx);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSCN);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL,(TSIFunction) IFunction,ctx);CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,J,J,(TSIJacobian)IJacobian,ctx);CHKERRQ(ierr);
  ierr = TSSetApplicationContext(ts,ctx);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetInitialGuess(X,ctx);CHKERRQ(ierr);

  ierr = VecDuplicate(X,&F_alg);CHKERRQ(ierr);
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes_alg);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes_alg,F_alg,AlgFunction,ctx);CHKERRQ(ierr);
  ierr = MatZeroEntries(J);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes_alg,J,J,AlgJacobian,ctx);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(snes_alg,"alg_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes_alg);CHKERRQ(ierr);
  /* Solve the algebraic equations */
  ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);

  ierr = SetSolution(ctx, ctx->t0, X); CHKERRQ(ierr);
  ierr = SetObservation(ctx, ctx->t0, X); CHKERRQ(ierr);


  /* Just to set up the Jacobian structure */
  ierr = VecDuplicate(X,&Xdot);CHKERRQ(ierr);
  ierr = IJacobian(ts,ctx->t0,X,Xdot,0.0,J,J,ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&Xdot);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve from on [t0,tdisturb] (steady state)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSetDuration(ts,1000000,ctx->tdisturb);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,ctx->t0,ctx->dt);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,SaveObservation);CHKERRQ(ierr);
  /* Solve (from t0 to tdisturb) */
  ierr = TSSolve(ts,X);CHKERRQ(ierr);

  /* set X at tdisturb as IC for the  optimization/estimation */
  /*ierr = VecDuplicate(X, &ctx->X0_disturb);CHKERRQ(ierr);*/
  /*ierr = VecCopy(X, ctx->X0_disturb);CHKERRQ(ierr);*/

  /* Continue integrating the DAE only if observations file is not specified */


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve from on [tdisturb, trestore] (disturbance part of the transient)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* Induce a load perturbation at t=tdisturb */
  for(i=0; i<3; i++) PD0[i] = PD0_disturb[i];
  
  printf("Generate data - initiated a bump in load: new PD0[0]=%g\n", PD0[0]);
  printf("Running with: tfinal=%.12f dt=%.12f data_dt=%.12f data_noise=%.12f prior_noise=%.12f load_disturb=%.12f\n",
	 ctx->tfinal, ctx->dt, ctx->data_dt,
	 ctx->data_noise, ctx->prior_noise, PD0[0]);

  /* Solve the algebraic equations  */
  ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);
  
  ierr = TSSetDuration(ts,100000,fmin(ctx->trestore,ctx->tfinal));CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,ctx->tdisturb,ctx->dt);CHKERRQ(ierr);
  /* Solve (from tdisturb to trestore) */
  ierr = TSSolve(ts,X);CHKERRQ(ierr);

  /* set X at trestore as IC for the  optimization/estimation */
  ierr = VecDuplicate(X, &ctx->X0_disturb);CHKERRQ(ierr);
  ierr = VecCopy(X, ctx->X0_disturb);CHKERRQ(ierr);

  
  if(0==strlen(ctx->loadObsFile)) {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Solve from on [trestore, tfinal] (post-disturbance transient)
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if(ctx->tfinal>=ctx->trestore+1.0e-8) {
      //restore  load at trestore
      for(i=0; i<3; i++) PD0[i] = PD0_ref[i];
      
      printf("In generate data: Restore load to PD0[0]=%g\n", PD0[0]);
      
      /* Solve the algebraic equations  */
      ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);
      
      ierr = TSSetDuration(ts,100000,ctx->tfinal);CHKERRQ(ierr);
      ierr = TSSetInitialTimeStep(ts,ctx->trestore,ctx->dt);CHKERRQ(ierr);
      /* Solve (from trestore to tfinal) */
      ierr = TSSolve(ts,X);CHKERRQ(ierr);
    } else {
      printf("Ignoring trestore since tfinal is less than it.\n");
    }
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Generate noise at level 'noise' percent
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = MatGetSize(ctx->obs, &m, &n);CHKERRQ(ierr);
    /* allocate std dev for data */
    ierr = PetscMalloc(m*sizeof(PetscReal),&ctx->data_stddev);CHKERRQ(ierr);
    ierr = MatDenseGetArray(ctx->obs,&mat);CHKERRQ(ierr);
    
    for(i=0;i<m;i++)
      ctx->data_stddev[i] = ctx->data_noise;
    
    /*  for(i=0; i<m; i++) {
	temp = 0.0;
	for(j=0; j<n; j++) {
	temp += mat[i*n+j]*mat[i*n+j];
	}
	ctx->data_stddev[i] = ctx->data_noise*sqrt(temp);
	}
    */
    
    for(i=0; i<m; i++) {
      for(j=0; j<n; j++) {
	mat[i*n+j] +=  ctx->data_stddev[i]*nrand();
      }
    }
    
    ierr = MatDenseRestoreArray(ctx->obs,&mat);CHKERRQ(ierr);
  } else {
    /* observations are in an external file */
    printf("Loading observations from %s.\n", ctx->loadObsFile);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,ctx->loadObsFile,FILE_MODE_READ, &obsView);CHKERRQ(ierr);
    ierr = MatLoad(ctx->obs, obsView);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&obsView);CHKERRQ(ierr);

    ierr = MatGetSize(ctx->obs, &m, &n);CHKERRQ(ierr);
    ierr = PetscMalloc(m*sizeof(PetscReal),&ctx->data_stddev);CHKERRQ(ierr);
    for(i=0;i<m;i++)
      ctx->data_stddev[i] = ctx->data_noise;
  }

  ierr = SNESDestroy(&snes_alg);CHKERRQ(ierr);
  ierr = VecDestroy(&F_alg);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------ */
/* Noemi: added this function, needs update!                          */
/* ------------------------------------------------------------------ */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionGradient"
/*
   FormFunction - Evaluates the function and corresponding gradient.

   Input Parameters:
   tao - the Tao context
   X   - the input vector
   ptr - optional user-defined context, as set by TaoSetObjectiveAndGradientRoutine()

   Output Parameters:
   f   - the newly evaluated function
   G   - the newly evaluated gradient
*/
PetscErrorCode FormFunctionGradient(Tao tao,Vec P,PetscReal *f,Vec G,void *ctx0)
{
  TS             ts;
  PetscErrorCode ierr;
  Userctx        *ctx = (Userctx*)ctx0;
  Vec            X, F_alg;
  SNES           snes_alg;
  PetscScalar    *x_ptr;
  Vec            lambda[1];
  //Vec            q;
  Vec            mu[1];
  PetscInt       steps,steps3;
  PetscReal      t,t2;
  Vec            Xdot;
  /* FD check */
  PetscReal      f1,f2,expo;
  Vec            Pvec_eps;
  PetscReal*     P_eps;
  PetscInt i;
  PetscBool fd;
  Vec Xdist_final;

  printf("aaa\n");

  ierr  = VecGetArray(P,&x_ptr);CHKERRQ(ierr);
  H[0] = x_ptr[0];
  H[1] = x_ptr[1];
  H[2] = x_ptr[2];
  //printf("FormFunctionGradient: x=[%.14f, %.14f, %.14f]\n",  x_ptr[0],  x_ptr[1], x_ptr[2]);
  //printf("FormFunctionGradient - PD0[0]=%g\n", PD0[0]);
  ierr  = VecRestoreArray(P,&x_ptr);CHKERRQ(ierr);

  if(ctx->t0 > ctx->tdisturb) {
    printf("t0 cannot be greater than tdisturb\n");
    PetscFunctionReturn(-1);
  }
  if( (ctx->tdisturb >= ctx->trestore-1.0e-8) || (ctx->tdisturb >= ctx->tfinal-1.0e-8) ) {
    printf("tdisturb should be less than trestore and tfinal\n");
    PetscFunctionReturn(-1);
  }

  ctx->misfit=0.0;
  ctx->stepnum = 0;

  ierr = VecZeroEntries(ctx->vec_q);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dmpgrid,&X);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSCN);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL,(TSIFunction) IFunction,ctx);CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,ctx->J,ctx->J,(TSIJacobian)IJacobian,ctx);CHKERRQ(ierr);
  ierr = TSSetApplicationContext(ts,ctx);CHKERRQ(ierr);

  /* Set initial conditions */
  ierr = VecCopy(ctx->X0_disturb, X);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve from on [tdisturb, trestore] (disturbance part of the transient)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* Induce a load perturbation at t=tdisturb */
  //!for(i=0; i<3; i++) PD0[i] = PD0_disturb[i];

  /* Induce a load perturbation at t=trestore*/
  for(i=0; i<3; i++) PD0[i] = PD0_ref[i];
  //!printf("In FormFunctionGradien: Induce a load perturbance to PD0[0]=%g\n", PD0[0]);

  /* Solve for algebraic variables with Xgen given by X0_disturb */
  ierr = VecDuplicate(X,&F_alg);CHKERRQ(ierr);
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes_alg);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes_alg,F_alg,AlgFunction,ctx);CHKERRQ(ierr);
  ierr = MatZeroEntries(ctx->J);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes_alg,ctx->J,ctx->J,AlgJacobian,ctx);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(snes_alg,"alg_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes_alg);CHKERRQ(ierr);
  /* Solve the algebraic equations */
  ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);

  /* Just to set up the Jacobian structure */
  ierr = VecDuplicate(X,&Xdot);CHKERRQ(ierr);
  //!  ierr = IJacobian(ts,ctx->tdisturb,X,Xdot,0.0,ctx->J,ctx->J,ctx);CHKERRQ(ierr);
  ierr = IJacobian(ts,ctx->trestore,X,Xdot,0.0,ctx->J,ctx->J,ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&Xdot);CHKERRQ(ierr);

  /* Save trajectory of solution so that TSAdjointSolve() may be used */
  ierr = TSSetSaveTrajectory(ts);CHKERRQ(ierr);

  /* Hook up the function evaluation */
  ierr = TSSetPostStep(ts,EvalMisfit);CHKERRQ(ierr);

  //!ierr = TSSetDuration(ts,10000,fmin(ctx->trestore,ctx->tfinal));CHKERRQ(ierr);
  ierr = TSSetDuration(ts,10000,ctx->tfinal);CHKERRQ(ierr);
  //!ierr = TSSetInitialTimeStep(ts,ctx->tdisturb,ctx->dt);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,ctx->trestore,ctx->dt);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  /* Solve the forward problem */
  //printf("Forward solve...\n");
  ierr = TSSolve(ts,X);CHKERRQ(ierr);

  ierr = VecDuplicate(X, &Xdist_final);CHKERRQ(ierr);
  ierr = VecCopy(X, Xdist_final);CHKERRQ(ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve from on [trestore, tfinal] (post-disturbance transient)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* if(ctx->tfinal>=ctx->trestore+1.0e-8) { */
  /*   //restore  load at trestore */
  /*   for(i=0; i<3; i++) PD0[i] = PD0_ref[i]; */
    
  /*   printf("In FormFunctionGradien: Restore load to PD0[0]=%g\n", PD0[0]); */
    
  /*   /\* Solve the algebraic equations  *\/ */
  /*   ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr); */
    
  /*   ierr = TSSetDuration(ts,100000,ctx->tfinal);CHKERRQ(ierr); */
  /*   ierr = TSSetInitialTimeStep(ts,ctx->trestore,ctx->dt);CHKERRQ(ierr); */
  /*   /\* Solve (from trestore to tfinal) *\/ */
  /*   ierr = TSSolve(ts,X);CHKERRQ(ierr); */
  /* } else { */
  /*   printf("Ignoring trestore since tfinal is less than it.\n"); */
  /* } */




  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Adjoint model starts here
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSGetTimeStepNumber(ts,&steps3);CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,NULL);CHKERRQ(ierr);

  ierr = MatCreateVecs(ctx->J,&lambda[0],NULL);CHKERRQ(ierr);

  /*   Set initial conditions for the adjoint integration */
  ierr = VecZeroEntries(lambda[0]);CHKERRQ(ierr);

  ierr = MatCreateVecs(ctx->Jacp,&mu[0],NULL);CHKERRQ(ierr);

  ierr = VecZeroEntries(mu[0]);CHKERRQ(ierr);

  /* Sets the initial value of the gradients of the cost w.r.t. x_0 and p */
  /*  Notes: the entries in these vectors must be correctly initialized */
  /* with the values lambda_i = df/dy|finaltime mu_i = df/dp|finaltime */
  ierr = TSSetCostGradients(ts,1,lambda,mu);CHKERRQ(ierr);

  /* Sets the function that computes the Jacobian of f w.r.t. p where x_t = f(x,y,p,t) */
  ierr = TSAdjointSetRHSJacobian(ts,ctx->Jacp,RHSJacobianP,ctx);CHKERRQ(ierr);

  /* Sets the routine for evaluating the integral term in the cost */
  /*ierr = TSSetCostIntegrand(ts,1,
			    (PetscErrorCode (*)(TS,PetscReal,Vec,Vec,void*))CostIntegrand,
			    (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DRDYFunction,
			    (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DRDPFunction,ctx);
  */
  ierr = TSSetCostIntegrand(ts,1,
			    NULL,
			    (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DRDYFunction,
			    (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DRDPFunction,ctx);
  CHKERRQ(ierr);

  t = ctx->tfinal;
  steps = (PetscInt)round(ctx->data_dt/ctx->dt);
  while(fabs(t-ctx->trestore)>1e-8)
  {
    ierr = TSGetTime(ts, &t2);CHKERRQ(ierr);

    /* Induce the perturbation in load accordingly corresponding to this time */
    if(t2-ctx->trestore>=-1e-8)
      for(i=0; i<3; i++) PD0[i] = PD0_ref[i];
    /* else if(t2-ctx->tdisturb>=0) */
    /*   for(i=0; i<3; i++) PD0[i] = PD0_disturb[i]; */
    else {printf("Panic: should not get here\n"); PetscFunctionReturn(-1);}

    /* Initial conditions for the adjoint */
    /* lambda += dr/dy */
    ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);
          
    ierr = AddDRDY(t2,X,&lambda[0],ctx);CHKERRQ(ierr);
    
    //printf("Manual adjoint backward integration steps=%d t=%g t2=%g \n", steps, t, t2);
    /* Sets # steps the adjoint solver should take backward in time*/
    ierr = TSAdjointSetSteps(ts,steps);CHKERRQ(ierr);

    /* Solves the discrete adjoint problem for an ODE/DAE */
    ierr = TSAdjointSolve(ts);CHKERRQ(ierr);

    t -= steps * ctx->dt;
  }

  //printf("mu-FunctionGradient after Adjoint (t=%g)\n",t);
  //ierr = VecView(mu[0],PETSC_VIEWER_STDOUT_SELF);
  //ierr = VecView(lambda[0],PETSC_VIEWER_STDOUT_SELF);

  /* return gradient */
  ierr = VecCopy(mu[0],G);CHKERRQ(ierr);
  ierr = AddRegGradient(ctx,P,G);

  //ierr = VecView(G,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);

  /* return fcn eval */
  *f  = ctx->misfit;
  EvalReg(ctx, P);
  *f += ctx->prior;
  //printf("objective=%.12f\n", *f);
  
  /* Finalize: destroy */
  ierr = VecDestroy(&lambda[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&mu[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&F_alg);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes_alg);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  //printf("Adjoint ends\n");

  fd=0;
  if(fd) {
    /* FD check */
    ierr =  FormFunction(tao,P,&f1,ctx); CHKERRQ(ierr);
    printf("cost=%.12f \n",f1);
    ierr = VecDuplicate(P, &Pvec_eps); CHKERRQ(ierr);

    for(i=0; i<3; i++) {
      for(expo=1e-2; expo>1e-8; expo/=3) {

	ierr = VecCopy(P, Pvec_eps); CHKERRQ(ierr);

	ierr = VecGetArray(Pvec_eps, &P_eps); CHKERRQ(ierr);

	P_eps[i] += expo;
	ierr = VecRestoreArray(Pvec_eps, &P_eps); CHKERRQ(ierr);

	//ierr = VecView(Pvec_eps,PETSC_VIEWER_STDOUT_SELF);

	ierr =  FormFunction(tao,Pvec_eps,&f2,ctx); CHKERRQ(ierr);
	printf("fd[%d]=%12.6e f1=%.7e f2=%.7e expo=%g\n", i+1, (f2-f1)/expo, f1, f2, expo);
      }
    }
    ierr = VecDestroy(&Pvec_eps); CHKERRQ(ierr); 
    /* ~end of FD */
  }
  //PetscFunctionReturn(-1);
  PetscFunctionReturn(0);
}
