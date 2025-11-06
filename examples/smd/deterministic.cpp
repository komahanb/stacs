#include "common.h"

#include <cstdio>
#include <vector>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsScalar mass = 2.5;
  TacsScalar damping = 0.2;
  TacsScalar stiffness = 5.0;
  TacsScalar parameters[3] = {mass, damping, stiffness};

  TacsScalar fvals[SMD_NUM_FUNCS];
  TacsScalar dfdx[SMD_NUM_FUNCS*SMD_NUM_DVARS];
  evaluate_smd_system(comm, parameters, true, 50.0, fvals, dfdx);

  printf("pe = %.17e, u = %.17e", TacsRealPart(fvals[0]), TacsRealPart(fvals[1]));
  printf("d{pe}dm = %.17e %.17e", TacsRealPart(dfdx[0]), TacsRealPart(dfdx[1]));
  printf("d{u}dm  = %.17e %.17e", TacsRealPart(dfdx[SMD_NUM_DVARS]), TacsRealPart(dfdx[SMD_NUM_DVARS+1]));

  const double dh = 1.0e-10;
  TacsScalar fhvals[SMD_NUM_FUNCS];
  TacsScalar dfdx_tmp[SMD_NUM_FUNCS*SMD_NUM_DVARS];
  TacsScalar perturbed_params[3];

  perturbed_params[0] = mass + dh;
  perturbed_params[1] = damping;
  perturbed_params[2] = stiffness;
  evaluate_smd_system(comm, perturbed_params, true, 50.0, fhvals, dfdx_tmp);
  printf("df1dm %.17e", TacsRealPart(fhvals[0]-fvals[0])/dh);
  printf("df2dm %.17e", TacsRealPart(fhvals[1]-fvals[1])/dh);

  perturbed_params[0] = mass;
  perturbed_params[1] = damping;
  perturbed_params[2] = stiffness + dh;
  evaluate_smd_system(comm, perturbed_params, true, 50.0, fhvals, dfdx_tmp);
  printf("df1dk %.17e", TacsRealPart(fhvals[0]-fvals[0])/dh);
  printf("df2dk %.17e", TacsRealPart(fhvals[1]-fvals[1])/dh);

  MPI_Finalize();
  return 0;
}
