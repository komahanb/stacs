#include "common.h"
#include "smd.h"

#include <cstdio>
#include <cstring>

#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "TACSIntegrator.h"

#include "TACSFunction.h"

#include "TACSPotentialEnergy.h"
#include "TACSDisplacement.h"
#include "TACSKSFunction.h"

#include "ParameterContainer.h"
#include "ParameterFactory.h"

#include "TACSStochasticElement.h"

#include "TACSStochasticFunction.h"
#include "TACSKSStochasticFunction.h"

#include "TACSStochasticFMeanFunction.h"
#include "TACSStochasticFFMeanFunction.h"

#include <vector>

void updateElement(TACSElement *elem, TacsScalar *vals, void *ctx){
  (void)ctx;
  SMD *smd = dynamic_cast<SMD*>(elem);
  if (smd != NULL) {
    smd->c = vals[0];
  } else {
    printf("Element mismatch while updating...");
  }
}

int main( int argc, char *argv[] ){

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank; 
  MPI_Comm_rank(comm, &rank); 

  FILE *fp = fopen("projection-smd-ks.dat", "w");

  for (int ii = 0; ii < 10; ii++){
    ParameterFactory *factory = new ParameterFactory();
    AbstractParameter *c = factory->createNormalParameter(0.2, 0.1, ii);

    ParameterContainer *pc = new ParameterContainer();
    pc->addParameter(c);
    pc->initialize();

    int nsterms = pc->getNumBasisTerms();

    TacsScalar mass = 2.5;
#ifdef TACS_USE_COMPLEX
    mass += TacsScalar(0.0, 1.0e-30);
#endif
    TacsScalar damping = 0.2;
    TacsScalar stiffness = 5.0;
    TACSElement *smd = new SMD(mass, damping, stiffness, 1.0, 0.0); 
    TACSStochasticElement *ssmd = new TACSStochasticElement(smd, pc, updateElement);

    // Assembler information to create TACS  
    int nelems = 1;
    int nnodes = 1;  
    int vars_per_node = nsterms;

    int conn[1] = {0};
    int ptr[2] = {0, 1};
    int eids[1] = {0};
    TacsScalar X[3] = {0.0, 0.0, 0.0};

    TACSCreator *creator = new TACSCreator(comm, vars_per_node);
    creator->incref();
    if (rank == 0){    
      creator->setGlobalConnectivity(nnodes, nelems, ptr, conn, eids);
      creator->setNodes(X);
    }
    TACSElement *elems[1] = { ssmd };
    creator->setElements(elems, nelems);

    TACSAssembler *tacs = creator->createTACS();
    tacs->incref();  
    creator->decref(); 

    TacsScalar dvs[SMD_NUM_DVARS] = { mass, stiffness };
    tacs->setDesignVars(dvs, SMD_NUM_DVARS);

    //---------------------------------------------------------------//  
    // Setup function evaluation within TACS
    //---------------------------------------------------------------//
  
    const int num_dvars = SMD_NUM_DVARS;
    const int num_funcs = 4;
    const int ks = 1;
    double ksweight = 10000.0;
      
    TACSFunction *pe, *disp;
    if (!ks){
      pe = new TACSPotentialEnergy(tacs);
      disp = new TACSDisplacement(tacs);
    } else {
      pe = new TACSKSFunction(tacs, TACS_POTENTIAL_ENERGY_FUNCTION, ksweight);
      disp = new TACSKSFunction(tacs, TACS_DISPLACEMENT_FUNCTION, ksweight);
    }

    TACSFunction *spe, *sdisp;
    TACSFunction *spe2, *sdisp2;
    if (!ks){
      spe = new TACSStochasticFMeanFunction(tacs, pe, pc, TACS_POTENTIAL_ENERGY_FUNCTION, FUNCTION_MEAN);
      sdisp = new TACSStochasticFMeanFunction(tacs, disp, pc, TACS_DISPLACEMENT_FUNCTION, FUNCTION_MEAN);

      spe2 = new TACSStochasticFFMeanFunction(tacs, pe, pc, TACS_POTENTIAL_ENERGY_FUNCTION, FUNCTION_VARIANCE);
      sdisp2 = new TACSStochasticFFMeanFunction(tacs, disp, pc, TACS_DISPLACEMENT_FUNCTION, FUNCTION_VARIANCE);
    } else {    
      spe = new TACSKSStochasticFunction(tacs, pe, pc, TACS_POTENTIAL_ENERGY_FUNCTION, FUNCTION_MEAN, ksweight);
      sdisp = new TACSKSStochasticFunction(tacs, disp, pc, TACS_DISPLACEMENT_FUNCTION, FUNCTION_MEAN, ksweight);

      spe2 = new TACSKSStochasticFunction(tacs, pe, pc, TACS_POTENTIAL_ENERGY_FUNCTION, FUNCTION_VARIANCE, ksweight);
      sdisp2 = new TACSKSStochasticFunction(tacs, disp, pc, TACS_DISPLACEMENT_FUNCTION, FUNCTION_VARIANCE, ksweight);
    }

    std::vector<TACSFunction*> funcs;
    funcs.push_back(spe);
    funcs.push_back(spe2);
    funcs.push_back(sdisp);
    funcs.push_back(sdisp2);
    for (auto *func : funcs){
      func->incref();
    }

    TacsScalar ftmp[num_funcs];
    memset(ftmp, 0, num_funcs*sizeof(TacsScalar));

    //-----------------------------------------------------------------//
    // Create the integrator class
    //-----------------------------------------------------------------//

    double tinit = 0.0;
    double tfinal = 10.0;
    int nsteps = 100;
    int time_order = 2;
    TACSIntegrator *bdf = new TACSBDFIntegrator(tacs, tinit, tfinal, nsteps, time_order);
    bdf->incref();
    bdf->setAbsTol(1e-12);
    bdf->setPrintLevel(0);
    bdf->setFunctions(funcs.data(), num_funcs, num_dvars);
    bdf->integrate();  
    bdf->evalFunctions(ftmp);
    bdf->integrateAdjoint();   

    TacsScalar pemean  = ftmp[0];
    TacsScalar pe2mean = ftmp[1];
    TacsScalar umean   = ftmp[2];
    TacsScalar u2mean  = ftmp[3];

    TacsScalar pevar = pe2mean - pemean*pemean;
    TacsScalar uvar  = u2mean - umean*umean;
    printf("Expectations : %.17e %.17e\n", RealPart(pemean), RealPart(umean));
    printf("Variance     : %.17e %.17e\n", RealPart(pevar), RealPart(uvar));

    std::vector<TacsScalar> grad(num_funcs * num_dvars, 0.0);
    bdf->getGradient(grad.data());

    TacsScalar pemeanderiv0 = grad[0*num_dvars + 0];
    TacsScalar pemeanderiv1 = grad[0*num_dvars + 1];
    TacsScalar pe2meanderiv0 = grad[1*num_dvars + 0];
    TacsScalar pe2meanderiv1 = grad[1*num_dvars + 1];
    TacsScalar umeanderiv0 = grad[2*num_dvars + 0];
    TacsScalar umeanderiv1 = grad[2*num_dvars + 1];
    TacsScalar u2meanderiv0 = grad[3*num_dvars + 0];
    TacsScalar u2meanderiv1 = grad[3*num_dvars + 1];

    TacsScalar pevar_deriv0 = pe2meanderiv0 - 2.0*pemean*pemeanderiv0;
    TacsScalar pevar_deriv1 = pe2meanderiv1 - 2.0*pemean*pemeanderiv1;
    TacsScalar uvar_deriv0 = u2meanderiv0 - 2.0*umean*umeanderiv0;
    TacsScalar uvar_deriv1 = u2meanderiv1 - 2.0*umean*umeanderiv1;

    double dh = 1.0e-30;
    printf("CS dE{ u  }/dx = %.17e %.17e %.17e\n", RealPart(umeanderiv0),
           ImagPart(umean)/dh,
           RealPart(umeanderiv0) - ImagPart(umean)/dh);
    printf("CS dE{ pe }/dx = %.17e %.17e %.17e\n", RealPart(pemeanderiv0),
           ImagPart(pemean)/dh,
           RealPart(pemeanderiv0) - ImagPart(pemean)/dh);
    printf("CS dV{ u  }/dx = %.17e %.17e %.17e\n", RealPart(uvar_deriv0),
           ImagPart(uvar)/dh,
           RealPart(uvar_deriv0) - ImagPart(uvar)/dh);
    printf("CS dV{ pe }/dx = %.17e %.17e %.17e\n", RealPart(pevar_deriv0),
           ImagPart(pevar)/dh,
           RealPart(pevar_deriv0) - ImagPart(pevar)/dh);

    fprintf(fp, "%d  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e\n",
            ii+1,
            RealPart(pemean),
            RealPart(pevar),
            RealPart(pemeanderiv0),
            RealPart(pevar_deriv0),
            ImagPart(pemean)/dh,
            ImagPart(pevar)/dh,    
            RealPart(pemeanderiv0) - ImagPart(pemean)/dh,
            RealPart(pevar_deriv0) - ImagPart(pevar)/dh);

    bdf->decref();
    for (auto *func : funcs){
      func->decref();
    }
    tacs->decref();
    ssmd->decref();
  }

  fclose(fp);

  MPI_Finalize();  
  return 0;
}
