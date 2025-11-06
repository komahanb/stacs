#include "TACSElement.h"

// Define some quantities of interest
static const int TACS_KINETIC_ENERGY_FUNCTION   = -1;
static const int TACS_POTENTIAL_ENERGY_FUNCTION = -2;
static const int TACS_DISPLACEMENT_FUNCTION     = -3;
static const int TACS_VELOCITY_FUNCTION         = -4;

class SMD : public TACSElement{  
 public:
  SMD(TacsScalar m, TacsScalar c, TacsScalar k, TacsScalar u0, TacsScalar udot0);
  ~SMD();

  /**
     Return the Initial conditions
  */
  void getInitConditions( TacsScalar v[], TacsScalar dv[],
                          TacsScalar ddv[], const TacsScalar X[] ) override;

  /**
     Compute the residual of the governing equations
  */
  void addResidual( double time,
                    TacsScalar res[],
                    const TacsScalar X[], const TacsScalar v[],
                    const TacsScalar dv[], const TacsScalar ddv[] ) override;

  /**
     Compute the Jacobian of the governing equations
  */
  void addJacobian( double time,
                    TacsScalar mat[],
                    double alpha, double beta, double gamma,
                    const TacsScalar X[], const TacsScalar v[],
                    const TacsScalar dv[], const TacsScalar ddv[] ) override;

  /**
     Evaluate a point-wise quantity of interest.
  */
  int evalPointQuantity( int elemIndex, int quantityType, double time,
                         int n, double pt[], const TacsScalar Xpts[],
                         const TacsScalar vars[], const TacsScalar dvars[],
                         const TacsScalar ddvars[], TacsScalar *quantity );

  void addPointQuantitySVSens( int elemIndex, int quantityType,
                               double time,
                               TacsScalar alpha,
                               TacsScalar beta,
                               TacsScalar gamma,
                               int n, double pt[],
                               const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfdq[],
                               TacsScalar dfdu[] );

  void addPointQuantityDVSens( int elemIndex, int quantityType,
                               double time,
                               TacsScalar scale,
                               int n, double pt[],
                               const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfdq[],
                               int dvLen,
                               TacsScalar dfdx[] );

  /**
     Get the element design variables values

     @param elemIndex The local element index
     @param dvLen The length of the design array
     @param dvs The design variable values
     @return The number of design variable numbers defined by the element
  */
  void getDesignVars( TacsScalar dvs[], int numDVs ) override{
    dvs[0] = this->m;
    dvs[1] = this->k;
  }

  /**
     Set the element design variables from the design vector

     @param elemIndex The local element index
     @param dvLen The length of the design array
     @param dvs The design variable values
     @return The number of design variable numbers defined by the element
  */
  void setDesignVars( const TacsScalar dvs[], int numDVs ) override{
    m = dvs[0];
    k = dvs[1];    
  }

  /**
     Get the lower and upper bounds for the design variable values

     @param elemIndex The local element index
     @param dvLen The length of the design array
     @param lowerBound The design variable lower bounds
     @param lowerBound The design variable upper bounds
     @return The number of design variable numbers defined by the element
  */
  void getDesignVarRange( TacsScalar lowerBound[],
                          TacsScalar upperBound[],
                          int numDVs ) override{
    // mass bounds
    lowerBound[0] = 1.0;
    upperBound[0] = 5.0;

    // stiffness bounds
    lowerBound[1] = 2.0;
    upperBound[1] = 10.0;
  }

  /**
     Add the derivative of the adjoint-residual product to the output vector

     This adds the contribution scaled by an input factor as follows:

     dvSens += scale*d(psi^{T}*(res))/dx

     By default the code is not implemented, but is not required so that
     analysis can be performed. Correct derivatives require a specific
     implementation.

     @param elemIndex The local element index
     @param time The simulation time
     @param scale The coefficient for the derivative result
     @param psi The element adjoint variables
     @param Xpts The element node locations
     @param vars The values of the element degrees of freedom
     @param dvars The first time derivative of the element DOF
     @param ddvars The second time derivative of the element DOF
     @param dvLen The length of the design variable vector
     @param dvSens The derivative vector
  */
  void addAdjResProduct( double time,
                         double scale,
                         TacsScalar dfdx[], int dvLen,
                         const TacsScalar psi[],
                         const TacsScalar Xpts[],
                         const TacsScalar vars[],
                         const TacsScalar dvars[],
                         const TacsScalar ddvars[] ) override;

  int numDisplacements() override{
    return 1;
  };
  
  int numNodes() override {
    return 1;
  }

  void setMass(TacsScalar m){
    // printf("updating mass [ %e -> %e ] \n", this->m, m);
    this->m = m;
  }

  void setStiffness(TacsScalar k){
    // printf("updating stiff [ %e -> %e ] \n", this->k, k);
    this->k = k;
  }

  void setDamping(TacsScalar c){
    // printf("updating damp [ %e -> %e ] \n", this->c, c);
    this->c = c;
  }

  void setInitPosition(TacsScalar u0){
    this->u0 = u0;
  }

  void setInitVelocity(TacsScalar udot0){
    this->udot0 = udot0;
  }

  // coefficients
  TacsScalar m, c, k;

  // initial conditions
  TacsScalar u0, udot0;
};
