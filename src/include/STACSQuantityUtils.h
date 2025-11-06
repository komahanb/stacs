#ifndef STACS_QUANTITY_UTILS_H
#define STACS_QUANTITY_UTILS_H

#include "TACSElement.h"

// Forward declaration of the SMD element to preserve existing behaviour
class SMD;

int STACSComputeQuantity( TACSElement *element,
                          int quantityType,
                          double time,
                          int n, double pt[],
                          const TacsScalar Xpts[],
                          const TacsScalar vars[],
                          const TacsScalar dvars[],
                          const TacsScalar ddvars[],
                          TacsScalar *quantity );

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
                             TacsScalar dfdu[] );

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
                             TacsScalar dfdx[] );

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
                              TacsScalar dfdX[] );

inline int STACSEvalPointQuantity( TACSElement *element,
                                   int elemIndex,
                                   int quantityType,
                                   double time,
                                   int n, double pt[],
                                   const TacsScalar Xpts[],
                                   const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   TacsScalar *quantity ){
  (void) elemIndex;
  return STACSComputeQuantity(element, quantityType, time,
                              n, pt, Xpts, vars, dvars, ddvars, quantity);
}

inline void STACSAddPointQuantitySVSens( TACSElement *element,
                                         int elemIndex,
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
  (void) elemIndex;
  STACSAddQuantitySVSens(element, quantityType, time,
                         alpha, beta, gamma,
                         n, pt, Xpts, vars, dvars, ddvars,
                         dfdq, dfdu);
}

inline void STACSAddPointQuantityDVSens( TACSElement *element,
                                         int elemIndex,
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
  (void) elemIndex;
  STACSAddQuantityDVSens(element, quantityType, time, scale,
                         n, pt, Xpts, vars, dvars, ddvars,
                         dfdq, dvLen, dfdx);
}

inline void STACSAddPointQuantityXptSens( TACSElement *element,
                                          int elemIndex,
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
  (void) elemIndex;
  STACSAddQuantityXptSens(element, quantityType, time, scale,
                          n, pt, Xpts, vars, dvars, ddvars,
                          dfdq, dfdX);
}

#endif
