#include "TACSElement.h"
class PPA : public TACSElement{  
 public:
  PPA( TacsScalar xcm, TacsScalar xf, 
       TacsScalar m  , TacsScalar If, 
       TacsScalar ch , TacsScalar ca,
       TacsScalar kh , TacsScalar ka );
  
  // Return the Initial conditions
  // -----------------------------
  void getInitConditions( int elemIndex, const TacsScalar X[],
                          TacsScalar v[], TacsScalar dv[], TacsScalar ddv[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( int elemIndex, double time,
                    const TacsScalar X[], const TacsScalar v[],
                    const TacsScalar dv[], const TacsScalar ddv[],
                    TacsScalar res[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( int elemIndex, double time,
                    TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                    const TacsScalar X[], const TacsScalar v[],
                    const TacsScalar dv[], const TacsScalar ddv[],
                    TacsScalar res[], TacsScalar mat[] );
  int getVarsPerNode(){
    return 2;
  };
  
  int getNumNodes() {
    return 1;
  }
  TacsScalar xcm, xf, m, If, ch, ca, kh, ka, s;
};
