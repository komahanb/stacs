#include "STACSQuantityUtils.h"
#include "smd.h"

#include <math.h>
#include <string.h>

static bool isSMD( TACSElement *element ){
  return dynamic_cast<SMD*>(element) != nullptr;
}

static SMD* asSMD( TACSElement *element ){
  return dynamic_cast<SMD*>(element);
}

int STACSComputeQuantity( TACSElement *element,
                          int quantityType,
                          double time,
                          int n, double pt[],
                          const TacsScalar Xpts[],
                          const TacsScalar vars[],
                          const TacsScalar dvars[],
                          const TacsScalar ddvars[],
                          TacsScalar *quantity ){
  if (!element || !quantity){
    return 0;
  }

  if (isSMD(element)){
    return asSMD(element)->evalPointQuantity(0, quantityType, time,
                                             n, pt, Xpts, vars, dvars, ddvars,
                                             quantity);
  }

  if (quantityType == TACS_KINETIC_ENERGY_FUNCTION ||
      quantityType == TACS_POTENTIAL_ENERGY_FUNCTION){
    TacsScalar Te = 0.0, Pe = 0.0;
    element->computeEnergies(time, &Te, &Pe, Xpts, vars, dvars);
    if (quantityType == TACS_KINETIC_ENERGY_FUNCTION){
      *quantity = Te;
    } else {
      *quantity = Pe;
    }
    return 1;
  }
  else if (quantityType == TACS_DISPLACEMENT_FUNCTION){
    int numVars = element->numVariables();
    TacsScalar accum = 0.0;
    for (int i = 0; i < numVars; i++){
      accum += vars[i]*vars[i];
    }
    *quantity = sqrt(accum);
    return 1;
  }
  else if (quantityType == TACS_VELOCITY_FUNCTION){
    int numVars = element->numVariables();
    TacsScalar accum = 0.0;
    for (int i = 0; i < numVars; i++){
      accum += dvars[i]*dvars[i];
    }
    *quantity = sqrt(accum);
    return 1;
  }

  // Default: quantity not supported
  *quantity = 0.0;
  return 0;
}

void STACSAddQuantitySVSens( TACSElement *element,
                             int quantityType,
                             double time,
                             double alpha, double beta, double gamma,
                             int n, double pt[],
                             const TacsScalar Xpts[],
                             const TacsScalar vars[],
                             const TacsScalar dvars[],
                             const TacsScalar ddvars[],
                             const TacsScalar dfdq[],
                             TacsScalar dfdu[] ){
  if (!element || !dfdu){
    return;
  }
  int numVars = element->numVariables();
  if (dfdu){
    memset(dfdu, 0, numVars*sizeof(TacsScalar));
  }

  if (isSMD(element)){
    asSMD(element)->addPointQuantitySVSens(0, quantityType, time,
                                           alpha, beta, gamma,
                                           n, pt, Xpts, vars, dvars, ddvars,
                                           dfdq, dfdu);
    return;
  }

  // For the L2 norms defined above, derivatives correspond to scaled vectors
  if (quantityType == TACS_DISPLACEMENT_FUNCTION){
    TacsScalar magnitude = 0.0;
    for (int i = 0; i < numVars; i++){
      magnitude += vars[i]*vars[i];
    }
    magnitude = sqrt(magnitude);
    if (magnitude > 1e-16){
      for (int i = 0; i < numVars; i++){
        dfdu[i] = alpha * vars[i] / magnitude;
      }
    }
    return;
  }

  if (quantityType == TACS_VELOCITY_FUNCTION){
    TacsScalar magnitude = 0.0;
    for (int i = 0; i < numVars; i++){
      magnitude += dvars[i]*dvars[i];
    }
    magnitude = sqrt(magnitude);
    if (magnitude > 1e-16){
      for (int i = 0; i < numVars; i++){
        dfdu[i] = beta * dvars[i] / magnitude;
      }
    }
    return;
  }
}

void STACSAddQuantityDVSens( TACSElement *element,
                             int quantityType,
                             double time,
                             double scale,
                             int n, double pt[],
                             const TacsScalar Xpts[],
                             const TacsScalar vars[],
                             const TacsScalar dvars[],
                             const TacsScalar ddvars[],
                             const TacsScalar dfdq[],
                             int dvLen,
                             TacsScalar dfdx[] ){
  if (!element || !dfdx){
    return;
  }

  if (isSMD(element)){
    asSMD(element)->addPointQuantityDVSens(0, quantityType, time,
                                           scale, n, pt,
                                           Xpts, vars, dvars, ddvars,
                                           dfdq, dvLen, dfdx);
    return;
  }
  // Default to zero for general elements
}

void STACSAddQuantityXptSens( TACSElement *element,
                              int quantityType,
                              double time,
                              double scale,
                              int n, double pt[],
                              const TacsScalar Xpts[],
                              const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[],
                              TacsScalar dfdX[] ){
  if (!element || !dfdX){
    return;
  }
  int numNodes = element->numNodes();
  memset(dfdX, 0, 3*numNodes*sizeof(TacsScalar));
  // Only SMD has explicit definitions (but no dependency on Xpts)
}
