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

static char help[] = "Driver for PCE.\n\n";
/* usage ./driver_pce.exe number_parameters parameters_file output_state */

#include "common-inpest.h"

PetscErrorCode PCEModelInit(Userctx *user);

#undef __FUNCT__
#define __FUNCT__ "DumpState"
PetscErrorCode DumpState(Userctx *user, Vec X)
{
  PetscErrorCode    ierr;
  Vec               Xgen,Xnet;
  PetscScalar       *xnet, *proj_vec;
  FILE              *file_out;
  PetscInt          obs_len;
  PetscInt          i;
  PetscFunctionBegin;

  ierr     = VecGetArray(user->proj,&proj_vec);CHKERRQ(ierr);
  ierr     = VecGetSize(user->proj, &obs_len);CHKERRQ(ierr);
  //idx      = 2*obs_len*(PetscInt) step_num;
  
  ierr     = DMCompositeGetLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  ierr     = DMCompositeScatter(user->dmpgrid,X,Xgen,Xnet);CHKERRQ(ierr);
  ierr     = VecGetArray(Xnet,&xnet);CHKERRQ(ierr);
  //ierr     = MatDenseGetArray(user->obs,&mat);CHKERRQ(ierr);
  
  file_out = fopen((char*)user->buf1, "a+");
  if(file_out==NULL) {
    printf("Couldn't open output file [%s] for writing.\n", (char*)user->buf1);
    PetscFunctionReturn(-1);
  }    
  
  for(i=0;i<obs_len; i++) {
    //user->misfit += 0.5*pow(xnet[2*((int)proj_vec[i])]   - mat[idx+2*i],  2); // obs[2*i, num_obs]
    //user->misfit += 0.5*pow(xnet[2*((int)proj_vec[i])+1] - mat[idx+2*i+1],2);  // obs[2*i+1, num_obs]
    
    fprintf(file_out, "%f %f ", xnet[2*((int)proj_vec[i])], xnet[2*((int)proj_vec[i])+1]);
  } 
  fclose(file_out);
  
  //ierr     = MatDenseRestoreArray(user->obs,&mat);CHKERRQ(ierr);
  ierr     = VecRestoreArray(Xnet,&xnet);CHKERRQ(ierr);
  ierr     = DMCompositeRestoreLocalVectors(user->dmpgrid,&Xgen,&Xnet);CHKERRQ(ierr);
  ierr     = VecRestoreArray(user->proj,&proj_vec);CHKERRQ(ierr);
  
  
  PetscFunctionReturn(0);
}

/* Evaluate the misfit at each observation time */
#undef __FUNCT__
#define __FUNCT__ "SaveForwardState"
PetscErrorCode SaveForwardState(TS ts)
{
  PetscErrorCode    ierr;
  Userctx           *user;
  Vec               X;
  PetscReal         t;

  PetscReal         step_num;

  PetscFunctionBegin;

  ierr     = TSGetApplicationContext(ts,&user);CHKERRQ(ierr);
  ierr     = TSGetTime(ts,&t);CHKERRQ(ierr);

  ierr     = TSGetSolution(ts,&X);CHKERRQ(ierr);

  //!t = t - user->tdisturb;
  t = t - user->trestore;
  step_num = round(t / user->data_dt);
  if(fabs(step_num - t/user->data_dt)<= 1e-6*user->dt) {
    ierr     = DumpState(user, X);CHKERRQ(ierr);
  }
  user->stepnum++;
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  FILE* file_param, *file_out;
  PetscErrorCode     ierr;
  PetscInt num_params, it, nread, num_simul;
  PetscReal val;
  Vec P;
  PetscReal* P_array, *proj_array;
  Userctx user;
  PetscInt i;

  ierr = PetscInitialize(&argc,&argv,"petscoptions",help);CHKERRQ(ierr);

  if(argc!=4) {
    printf("Invalid number of parameters. Usage: %s number_parameters parameters_file output_state\n", argv[0]);
    ierr = PetscFinalize();
    return (-1);
  }

  num_params = atoi(argv[1]);

  file_param = fopen(argv[2], "r+");
  if(NULL==file_param) {
    printf("Couldn't open file [%s] for reading.\n", argv[2]);
    ierr = PetscFinalize();
    return (-1);
  }

  user.buf1=(char*)malloc(sizeof(char)*strlen(argv[3]));
  strcpy(user.buf1,argv[3]);

  /* Model setup */
  ierr = ModelSetup(&user);CHKERRQ(ierr);

  /* Model initialization: simulate the systems to ensure steady state (1 second)
     and save the state to be used as IC for the forward simulations with different
     parameters (perturbed inertia). These simulations will be from t=1 to t=5.
  */
  ierr = PCEModelInit(&user);CHKERRQ(ierr);

  /* Parameter initialization */
  /* hard code the data projection/observable map - for now assume observations at all buses */
  ierr = VecCreateSeq(PETSC_COMM_WORLD,nbus,&user.proj);CHKERRQ(ierr);
  ierr = VecGetArray(user.proj,&proj_array);CHKERRQ(ierr);
  for(i=0; i<nbus; i++) {
    proj_array[i]=i;
  }
  ierr = VecRestoreArray(user.proj,&proj_array);CHKERRQ(ierr);

  /* create the parameters vector, load it from file and simulate forward for each */
  ierr = VecCreateSeq(PETSC_COMM_WORLD,num_params,&P);CHKERRQ(ierr);

  num_simul=0;
  while(1) {
    num_simul++;

    ierr = VecGetArray(P,&P_array);CHKERRQ(ierr);

    for(it=0; it<num_params; it++) {
      nread = fscanf(file_param, "%lf", &val); 
      
      if(nread==0) {
	printf("Error while reading the file due to a non-numeric value.\n");
	ierr = PetscFinalize();CHKERRQ(ierr);
	return (-1);
      } else if(nread<0) break;

      P_array[it]=val;
    }
    ierr = VecRestoreArray(P,&P_array);CHKERRQ(ierr);

    if(feof(file_param)) break;
    printf("Processing parameter vector [%3d]...", num_simul);  

    /* Solve the DAE and save the state */
    ierr = ForwardSolve(P,&user, SaveForwardState);CHKERRQ(ierr);

    //end with an end of line
    file_out = fopen((char*)user.buf1, "a+");
    if(file_out==NULL) {
      printf("Couldn't open output file [%s] for writing.\n", (char*)user.buf1);
      PetscFunctionReturn(-1);
    }
    fprintf(file_out, "\n");
    fclose(file_out);

    printf("done.\n");
  }  

  fclose(file_param);
  free(user.buf1);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCEModelInit"
PetscErrorCode PCEModelInit(Userctx *ctx)
{
  TS             ts;
  SNES           snes_alg;  
  Vec            X, F_alg, Xdot;
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;

  /* simulate with default inertia */
  H[0] = H0[0]; H[1] = H0[1]; H[2] = H0[2];

  ctx->stepnum = 0;

  //use the reference PD0 from t0 to t_disturb to ensure steady state
  for(i=0; i<3; i++) PD0[i] = PD0_ref[i];

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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetInitialGuess(X,ctx);CHKERRQ(ierr);

  ierr = VecDuplicate(X,&F_alg);CHKERRQ(ierr);
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes_alg);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes_alg,F_alg,AlgFunction,ctx);CHKERRQ(ierr);
  ierr = MatZeroEntries(ctx->J);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes_alg,ctx->J,ctx->J,AlgJacobian,ctx);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(snes_alg,"alg_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes_alg);CHKERRQ(ierr);
  /* Solve the algebraic equations */
  ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);
  ierr = VecDestroy(&F_alg);CHKERRQ(ierr);

  /* Just to set up the Jacobian structure */
  ierr = VecDuplicate(X,&Xdot);CHKERRQ(ierr);
  ierr = IJacobian(ts,0.0,X,Xdot,0.0,ctx->J,ctx->J,ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&Xdot);CHKERRQ(ierr);

  ctx->stepnum++;

  ierr = TSSetDuration(ts,10000,ctx->tdisturb);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,ctx->t0,ctx->dt);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  /*ierr = TSSetPostStep(ts,SaveObservation);CHKERRQ(ierr);*/
  ierr = TSSolve(ts,X);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve from on [tdisturb, trestore] (disturbance part of the transient)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* Induce a load perturbation at t=tdisturb */
  for(i=0; i<3; i++) PD0[i] = PD0_disturb[i];
  
  /* Solve the algebraic equations  */
  ierr = SNESSolve(snes_alg,NULL,X);CHKERRQ(ierr);
  
  ierr = TSSetDuration(ts,100000,fmin(ctx->trestore,ctx->tfinal));CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,ctx->tdisturb,ctx->dt);CHKERRQ(ierr);
  /* Solve (from tdisturb to trestore) */
  ierr = TSSolve(ts,X);CHKERRQ(ierr);

  // set X at tdisturb as IC for the  optimization/estimation
  ierr = VecDuplicate(X, &ctx->X0_disturb);CHKERRQ(ierr);
  ierr = VecCopy(X, ctx->X0_disturb);CHKERRQ(ierr);

  ierr = SNESDestroy(&snes_alg);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
