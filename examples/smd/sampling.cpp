#include"smd.h"

#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "TACSIntegrator.h"

#include "TACSFunction.h"
#include "TACSPotentialEnergy.h"
#include "TACSDisplacement.h"
#include "TACSKSFunction.h"

#include "ParameterContainer.h"
#include "ParameterFactory.h"

void deterministic_solve( MPI_Comm comm,
                          TacsScalar *p,
                          TacsScalar *fvals,
                          TacsScalar **dfdxvals ){
  int rank; 
  MPI_Comm_rank(comm, &rank); 

  TACSElement *smd = new SMD(p[0], p[1], p[2]); 
  smd->incref();

  // Assembler information to create TACS  
  int nelems = 1;
  int nnodes = 1;  
  int vars_per_node = 1;

  // Array of elements
  TACSElement **elems = new TACSElement*[nelems];
  elems[0] = smd;

  // Node points array
  TacsScalar *X = new TacsScalar[3*nnodes];
  memset(X, 0, 3*nnodes*sizeof(TacsScalar));

  // Connectivity array
  int *conn = new int[1];
  conn[0] = 0;

  // Connectivity pointer array
  int *ptr = new int[2];
  ptr[0] = 0;
  ptr[1] = 1;

  // Element Ids array
  int *eids = new int[nelems];
  for (int i = 0; i < nelems; i++){
    eids[i] = i;
  }

  // Creator object for TACS
  TACSCreator *creator = new TACSCreator(comm, vars_per_node);
  creator->incref();
  if (rank == 0){    
    creator->setGlobalConnectivity(nnodes, nelems, ptr, conn, eids);
    creator->setNodes(X);
  }
  creator->setElements(nelems, elems);

  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();  
  creator->decref();

  //---------------------------------------------------------------//  
  // Setup function evaluation within TACS
  //---------------------------------------------------------------//
  
  const int num_funcs = 2;
  TACSFunction *pe, *disp;
  int ks = 1;
  if (!ks) {
    pe = new TACSPotentialEnergy(tacs);
    disp = new TACSDisplacement(tacs);
  } else {
    double ksweight = 10000.0;
    pe = new TACSKSFunction(tacs, TACS_POTENTIAL_ENERGY_FUNCTION, ksweight);
    disp = new TACSKSFunction(tacs, TACS_DISPLACEMENT_FUNCTION, ksweight);
  }
  pe->incref();
  disp->incref();
  
  TACSFunction **funcs = new TACSFunction*[num_funcs];
  funcs[0] = pe;
  funcs[1] = disp;

  TACSBVec *dfdx1 = tacs->createDesignVec();  dfdx1->incref();
  TACSBVec *dfdx2 = tacs->createDesignVec();  dfdx2->incref();
  
  //---------------------------------------------------------------//
  // Create the integrator class
  //---------------------------------------------------------------//

  double tinit = 0.0;
  double tfinal = 10.0;
  int nsteps = 100;
  int time_order = 2;     
  TACSIntegrator *bdf = new TACSBDFIntegrator(tacs, tinit, tfinal, nsteps, time_order);
  bdf->incref();
  bdf->setAbsTol(1e-12);
  bdf->setPrintLevel(0);
  bdf->setFunctions(num_funcs, funcs);
  bdf->integrate();  
  bdf->evalFunctions(fvals);

  bdf->integrateAdjoint();
  bdf->getGradient(0, &dfdx1);
  bdf->getGradient(1, &dfdx2);

  TacsScalar *dfdx1vals;
  TacsScalar *dfdx2vals;
  dfdx1->getArray(&dfdx1vals);
  dfdx2->getArray(&dfdx2vals);

  const int num_dvars = 2;
  dfdxvals[0][0] = dfdx1vals[0];
  dfdxvals[0][1] = dfdx1vals[1];

  dfdxvals[1][0] = dfdx2vals[0];
  dfdxvals[1][1] = dfdx2vals[1];
  
  // clear allocated heap
  delete [] X;
  delete [] ptr;
  delete [] eids;
  delete [] conn;
  delete [] elems;
  delete [] funcs;

  pe->decref();
  disp->decref();  
  dfdx1->decref();
  dfdx2->decref();  
  smd->decref();
}


void sampling_solve(MPI_Comm comm,
                    int num_quad_pts,
                    TacsScalar mass, TacsScalar stiffness, 
                    TacsScalar *fmvals,
                    TacsScalar *fvvals,
                    TacsScalar *fmeanderiv=NULL, 
                    TacsScalar *fvarderiv=NULL){

  const int num_funcs = 2;
  const int num_dvars = 2;

  int pnqpts[1] = {num_quad_pts};
  ParameterFactory *factory = new ParameterFactory();
  AbstractParameter *c = factory->createNormalParameter(0.2, 0.1, 4);  

  ParameterContainer *pc = new ParameterContainer();
  pc->addParameter(c);
  pc->initializeQuadrature(pnqpts);

  const int nqpoints = pc->getNumQuadraturePoints();
  const int nvars = pc->getNumParameters();

  TacsScalar **f = new TacsScalar*[nqpoints];  
  for (int i = 0; i < nqpoints; i++){
    f[i] = new TacsScalar[num_funcs];  
  }

  TacsScalar ***dfdx = new TacsScalar**[nqpoints];
  for (int i = 0; i < nqpoints; i++){
    dfdx[i] = new TacsScalar*[num_funcs];
    for (int j = 0; j < num_funcs; j++){
      dfdx[i][j] = new TacsScalar[num_dvars];  
    }
  }
  
  TacsScalar *zq = new TacsScalar[nvars];
  TacsScalar *yq = new TacsScalar[nvars];
  TacsScalar wq;

  for (int q = 0; q < nqpoints; q++){
    wq = pc->quadrature(q, zq, yq);
    printf("deterministic solve %d at c = %.17e\n", q, RealPart(yq[0]));
    TacsScalar damping = yq[0];
    TacsScalar params[3] = {mass, damping, stiffness};
    deterministic_solve(comm, params, f[q], dfdx[q]);
    printf("\t disp = %.17e energy = %.17e \n", RealPart(f[q][0]), RealPart(f[q][1]));
  }

  //-------------------------------------------------------------------//
  // Post processing to find statistical moments (mean and variance)
  //-------------------------------------------------------------------//

  // Function values
  TacsScalar *fmean = new TacsScalar[num_funcs];  
  memset(fmean, 0, num_funcs*sizeof(TacsScalar));
  TacsScalar *fvar = new TacsScalar[num_funcs];  
  memset(fvar, 0, num_funcs*sizeof(TacsScalar));  

  // E[F*F] used to compute variance
  TacsScalar *f2mean = new TacsScalar[num_funcs];  
  memset(f2mean, 0, num_funcs*sizeof(TacsScalar));

  // Derivative values
  TacsScalar **dfdxmean = new TacsScalar*[num_funcs];
  for (int j = 0; j < num_funcs; j++){
    dfdxmean[j] = new TacsScalar[num_dvars];  
    memset(dfdxmean[j], 0, num_dvars*sizeof(TacsScalar));
  }
  TacsScalar **dfdxvar = new TacsScalar*[num_funcs];
  for (int j = 0; j < num_funcs; j++){
    dfdxvar[j] = new TacsScalar[num_dvars];  
    memset(dfdxvar[j], 0, num_dvars*sizeof(TacsScalar));
  }

  // Used for variance derivative
  TacsScalar **E2ffprime = new TacsScalar*[num_funcs];
  for (int j = 0; j < num_funcs; j++){
    E2ffprime[j] = new TacsScalar[num_dvars];  
    memset(E2ffprime[j], 0, num_dvars*sizeof(TacsScalar));
  }

  // Compute mean of function and derivatives
  for (int q = 0; q < nqpoints; q++){
    wq = pc->quadrature(q, zq, yq);
    for (int i = 0; i < num_funcs; i++){
      fmean[i] += wq*f[q][i]; // E[F]
      f2mean[i] += wq*f[q][i]*f[q][i]; // E[F*F]
      for (int j = 0; j < num_dvars; j++){
        dfdxmean[i][j] += wq*dfdx[q][i][j]; // E[dfdx]
        E2ffprime[i][j] += wq*2.0*f[q][i]*dfdx[q][i][j]; // E[2*F*dfdx]
      }
    }
  }

  // Print output expectation
  for (int i = 0; i < num_funcs; i++){
    printf("E[f%d] = %.17e \n", i, RealPart(fmean[i]));
    fmvals[i] = fmean[i];
  }  
  int idx = 0;
  for (int i = 0; i < num_funcs; i++){
    printf("E[df%ddx] = ", i);
    for (int j = 0; j < num_dvars; j++){
      printf("%.17e ", RealPart(dfdxmean[i][j]));
      if (fmeanderiv){ fmeanderiv[idx] = dfdxmean[i][j]; }
      idx++;
    }
    printf("\n");
  }
    
  // Compute variance of function
  for (int i = 0; i < num_funcs; i++){      
    fvar[i] = f2mean[i] - fmean[i]*fmean[i];
  }    
  // Compute variance derivatives
  for (int i = 0; i < num_funcs; i++){
    for (int j = 0; j < num_dvars; j++){
      dfdxvar[i][j] = E2ffprime[i][j] - 2.0*fmean[i]*dfdxmean[i][j];
    }
  }

  // Print output variance and variance derivative
  for (int i = 0; i < num_funcs; i++){
    printf("V[f%d] = %.17e \n", i, RealPart(fvar[i]));
    fvvals[i] = fvar[i];
  }  
  idx = 0;
  for (int i = 0; i < num_funcs; i++){
    printf("V[df%ddx] = ", i);
    for (int j = 0; j < num_dvars; j++){
      printf("%.17e ", RealPart(dfdxvar[i][j]));
      if (fvarderiv){ fvarderiv[idx] = dfdxvar[i][j]; }
      idx++;
    }
    printf("\n");
  }

  delete [] zq;
  delete [] yq;
   
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank; 
  MPI_Comm_rank(comm, &rank); 

  FILE *fp = fopen("sampling-smd-ks.dat", "w");
  for (int ii = 1; ii < 10; ii++){

  int num_quad_pts = ii;
  TacsScalar *fvals;           
  TacsScalar **dfdxvals;

  // Function values
  int num_funcs = 2; 
  int num_dvs = 2;

  TacsScalar *fmean = new TacsScalar[num_funcs];  
  memset(fmean, 0, num_funcs*sizeof(TacsScalar));
  TacsScalar *fvar = new TacsScalar[num_funcs];  
  memset(fvar, 0, num_funcs*sizeof(TacsScalar));  

  TacsScalar *fmeanderiv = new TacsScalar[num_funcs*num_dvs];  
  memset(fmeanderiv, 0, num_dvs*num_funcs*sizeof(TacsScalar));
  TacsScalar *fvarderiv = new TacsScalar[num_funcs*num_dvs];
  memset(fvarderiv, 0, num_dvs*num_funcs*sizeof(TacsScalar));

  double dh = 1.0e-30;
  TacsScalar mass = 2.5 + 1.0e-30j;
  TacsScalar stiffness = 5.0;

  // // FD Derivative of E[F] wrt mass
  // TacsScalar *fmeanmtmp = new TacsScalar[num_funcs];  
  // memset(fmeanmtmp, 0, num_funcs*sizeof(TacsScalar));
  // TacsScalar *fvarmtmp = new TacsScalar[num_funcs];  
  // memset(fvarmtmp, 0, num_funcs*sizeof(TacsScalar));  
  // sampling_solve(comm, mass + dh, stiffness, fmeanmtmp, fvarmtmp);

  // // FD Derivative of E[F] wrt stiffness
  // TacsScalar *fmeanktmp = new TacsScalar[num_funcs];  
  // memset(fmeanktmp, 0, num_funcs*sizeof(TacsScalar));
  // TacsScalar *fvarktmp = new TacsScalar[num_funcs];  
  // memset(fvarktmp, 0, num_funcs*sizeof(TacsScalar));  
  // sampling_solve(comm, mass, stiffness + dh, fmeanktmp, fvarktmp);
  
  // Baseline solution
  sampling_solve(comm, num_quad_pts, mass, stiffness, fmean, fvar, fmeanderiv, fvarderiv);
  TacsScalar pemean, pe2mean, pevar, pemeanderiv, pe2meanderiv;
  pemean = fmean[0];
  pevar = fvar[0];

  pemeanderiv = fmeanderiv[0];
  pe2meanderiv = fvarderiv[0];


  printf("Derivative of Expectation\n");
  int ctr = 0;
  for (int i = 0; i < num_funcs; i++){
    for (int j = 0; j < num_dvs; j++){
      printf("%d fd = %.17e actual = %.17e error = %.17e \n", i, 
             ImagPart(fmean[j])/dh, 
             RealPart(fmeanderiv[i+j*num_dvs]), 
             ImagPart(fmean[j])/dh - RealPart(fmeanderiv[i+j*num_dvs])
             );
      ctr++;
    }
  }

  printf("Derivative of Variance\n");
  ctr = 0;
  for (int i = 0; i < num_funcs; i++){
    for (int j = 0; j < num_dvs; j++){
      printf("%d fd = %.17e actual = %.17e error = %.17e \n", i, 
             ImagPart(fvar[j])/dh, 
             RealPart(fvarderiv[i+j*num_dvs]), 
             ImagPart(fvar[j])/dh - RealPart(fvarderiv[i+j*num_dvs])
             );
      ctr++;
    }
  }


  fprintf(fp, "%d  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e\n", 
          ii,
          RealPart(pemean),
          RealPart(pevar),
          RealPart(pemeanderiv),
          RealPart(pe2meanderiv),
          ImagPart(pemean)/dh,
          ImagPart(pevar)/dh,    
          RealPart(pemeanderiv) - ImagPart(pemean)/dh,
          RealPart(pe2meanderiv) - ImagPart(pevar)/dh);
    }


    fclose(fp);
  // printf("Derivative of Expectation\n");
  // int ctr = 0;
  // for (int i = 0; i < num_funcs; i++){
  //   for (int j = 0; j < num_dvs; j++){
  //     if (i == 0){
  //       printf("%d fd = %.17e actual = %.17e error = %.17e \n", i, 
  //              RealPart((fmeanmtmp[j] - fmean[j])/dh), 
  //              RealPart(fmeanderiv[i+j*num_dvs]), 
  //              RealPart((fmeanmtmp[j] - fmean[j])/dh - fmeanderiv[i+j*num_dvs]));
  //     } else { 
  //       printf("%d fd = %.17e actual = %.17e error = %.17e \n", i, 
  //              RealPart((fmeanktmp[j] - fmean[j])/dh), 
  //              RealPart(fmeanderiv[i+j*num_dvs]), 
  //              RealPart((fmeanktmp[j] - fmean[j])/dh - fmeanderiv[i+j*num_dvs]));
  //     } 
  //     ctr++;
  //   }
  // }
  /*  
  printf("Derivative of Variance\n");
  ctr = 0;
  for (int i = 0; i < num_funcs; i++){
    for (int j = 0; j < num_dvs; j++){
      if (i == 0){
        printf("%d fd = %.17e actual = %.17e error = %.17e \n", i, 
               RealPart((fvarmtmp[j] - fvar[j])/dh), 
               RealPart(fvarderiv[i+j*num_dvs]), 
               RealPart((fvarmtmp[j] - fvar[j])/dh - fvarderiv[i+j*num_dvs]));
      } else { 
        printf("%d fd = %.17e actual = %.17e error = %.17e \n", i, 
               RealPart((fvarktmp[j] - fvar[j])/dh), 
               RealPart(fvarderiv[i+j*num_dvs]), 
               RealPart((fvarktmp[j] - fvar[j])/dh - fvarderiv[i+j*num_dvs]));
      } 
      ctr++;
    }
  }
  */
  MPI_Finalize();  
  return 0;
}
