
#include "common.h"

#include "ParameterContainer.h"
#include "ParameterFactory.h"

#include <cstdio>
#include <vector>

void sampling_solve(MPI_Comm comm,
                    int num_quad_pts,
                    TacsScalar mass, TacsScalar stiffness,
                    TacsScalar *fmvals,
                    TacsScalar *fvvals,
                    TacsScalar *fmeanderiv = NULL,
                    TacsScalar *fvarderiv = NULL){
  const int num_funcs = SMD_NUM_FUNCS;
  const int num_dvars = SMD_NUM_DVARS;

  int pnqpts[1] = {num_quad_pts};
  ParameterFactory *factory = new ParameterFactory();
  AbstractParameter *c = factory->createNormalParameter(0.2, 0.1, 4);

  ParameterContainer *pc = new ParameterContainer();
  pc->addParameter(c);
  pc->initializeQuadrature(pnqpts);

  const int nqpoints = pc->getNumQuadraturePoints();
  const int nvars = pc->getNumParameters();

  std::vector<TacsScalar> fmean(num_funcs, 0.0);
  std::vector<TacsScalar> f2mean(num_funcs, 0.0);
  std::vector<TacsScalar> fvar(num_funcs, 0.0);

  std::vector<TacsScalar> dfdxmean(num_funcs*num_dvars, 0.0);
  std::vector<TacsScalar> dfdxvar(num_funcs*num_dvars, 0.0);
  std::vector<TacsScalar> E2ffprime(num_funcs*num_dvars, 0.0);

  std::vector<TacsScalar> fvals_local(num_funcs, 0.0);
  std::vector<TacsScalar> dfdx_local(num_funcs*num_dvars, 0.0);

  std::vector<TacsScalar> zq(nvars);
  std::vector<TacsScalar> yq(nvars);

  for (int q = 0; q < nqpoints; q++){
    TacsScalar wq = pc->quadrature(q, zq.data(), yq.data());
    printf("deterministic solve %d at c = %.17e\n", q, RealPart(yq[0]));
    TacsScalar damping = yq[0];
    TacsScalar params[3] = {mass, damping, stiffness};
    evaluate_smd_system(comm, params, true, 10000.0, fvals_local.data(), dfdx_local.data());
    printf("\t disp = %.17e energy = %.17e\n", RealPart(fvals_local[0]), RealPart(fvals_local[1]));

    for (int i = 0; i < num_funcs; i++){
      fmean[i] += wq*fvals_local[i];
      f2mean[i] += wq*fvals_local[i]*fvals_local[i];
      for (int j = 0; j < num_dvars; j++){
        int idx = i*num_dvars + j;
        dfdxmean[idx] += wq*dfdx_local[idx];
        E2ffprime[idx] += wq*2.0*fvals_local[i]*dfdx_local[idx];
      }
    }
  }

  for (int i = 0; i < num_funcs; i++){
    fvar[i] = f2mean[i] - fmean[i]*fmean[i];
    if (fvvals){ fvvals[i] = fvar[i]; }
    if (fmvals){ fmvals[i] = fmean[i]; }
  }

  for (int i = 0; i < num_funcs; i++){
    printf("E[f%d] = %.17e\n", i, RealPart(fmean[i]));
  }
  for (int i = 0; i < num_funcs; i++){
    printf("V[f%d] = %.17e\n", i, RealPart(fvar[i]));
  }

  for (int i = 0; i < num_funcs; i++){
    for (int j = 0; j < num_dvars; j++){
      int idx = i*num_dvars + j;
      dfdxvar[idx] = E2ffprime[idx] - 2.0*fmean[i]*dfdxmean[idx];
    }
  }

  if (fmeanderiv){
    for (int i = 0; i < num_funcs; i++){
      for (int j = 0; j < num_dvars; j++){
        int idx = i*num_dvars + j;
        fmeanderiv[idx] = dfdxmean[idx];
      }
    }
  }
  if (fvarderiv){
    for (int i = 0; i < num_funcs; i++){
      for (int j = 0; j < num_dvars; j++){
        int idx = i*num_dvars + j;
        fvarderiv[idx] = dfdxvar[idx];
      }
    }
  }

  for (int i = 0; i < num_funcs; i++){
    printf("E[df%ddx] = ", i);
    for (int j = 0; j < num_dvars; j++){
      int idx = i*num_dvars + j;
      printf("%.17e ", RealPart(dfdxmean[idx]));
    }
    printf("\n");
  }

  for (int i = 0; i < num_funcs; i++){
    printf("V[df%ddx] = ", i);
    for (int j = 0; j < num_dvars; j++){
      int idx = i*num_dvars + j;
      printf("%.17e ", RealPart(dfdxvar[idx]));
    }
    printf("\n");
  }
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank; 
  MPI_Comm_rank(comm, &rank); 

  FILE *fp = fopen("sampling-smd-ks.dat", "w");
  for (int ii = 1; ii < 10; ii++){
    int num_quad_pts = ii;

    std::vector<TacsScalar> fmean(SMD_NUM_FUNCS, 0.0);
    std::vector<TacsScalar> fvar(SMD_NUM_FUNCS, 0.0);
    std::vector<TacsScalar> fmeanderiv(SMD_NUM_FUNCS*SMD_NUM_DVARS, 0.0);
    std::vector<TacsScalar> fvarderiv(SMD_NUM_FUNCS*SMD_NUM_DVARS, 0.0);

    double dh = 1.0e-30;
    TacsScalar mass = 2.5;
#ifdef TACS_USE_COMPLEX
    mass += TacsScalar(0.0, dh);
#endif
    TacsScalar stiffness = 5.0;

    sampling_solve(comm, num_quad_pts, mass, stiffness,
                   fmean.data(), fvar.data(),
                   fmeanderiv.data(), fvarderiv.data());

    TacsScalar pemean = fmean[0];
    TacsScalar pevar = fvar[0];
    TacsScalar pemeanderiv = fmeanderiv[0];
    TacsScalar pe2meanderiv = fvarderiv[0];

    printf("Derivative of Expectation\n");
    for (int i = 0; i < SMD_NUM_FUNCS; i++){
      for (int j = 0; j < SMD_NUM_DVARS; j++){
        int idx = i*SMD_NUM_DVARS + j;
        printf("%d fd = %.17e actual = %.17e error = %.17e\n", i,
               ImagPart(fmean[j])/dh,
               RealPart(fmeanderiv[idx]),
               ImagPart(fmean[j])/dh - RealPart(fmeanderiv[idx]));
      }
    }

    printf("Derivative of Variance\n");
    for (int i = 0; i < SMD_NUM_FUNCS; i++){
      for (int j = 0; j < SMD_NUM_DVARS; j++){
        int idx = i*SMD_NUM_DVARS + j;
        printf("%d fd = %.17e actual = %.17e error = %.17e\n", i,
               ImagPart(fvar[j])/dh,
               RealPart(fvarderiv[idx]),
               ImagPart(fvar[j])/dh - RealPart(fvarderiv[idx]));
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

  MPI_Finalize();
  return 0;
}
