#ifndef TACS_STOCHASTIC_ELEMENT
#define TACS_STOCHASTIC_ELEMENT

#include "TACSElement.h"
#include "ParameterContainer.h"
#include "Python.h"

class TACSStochasticElement : public TACSElement {
 public:
  TACSStochasticElement( TACSElement *_delem,
                         ParameterContainer *_pc,
                         void (*_update)(TACSElement*, TacsScalar*, void*) );
  ~TACSStochasticElement();

  void setPythonCallback(PyObject *cbptr){
    this->pyptr = cbptr;
  }

  // TACS Element member functions
  // -----------------------------
  int numDisplacements() override;
  int numNodes() override;
  int numVariables() override;

  // Get the element basis
  //-----------------------
  void getMultiplierIndex( int *multiplier ) override{
    delem->getMultiplierIndex(multiplier);
  }

  // Return the Initial conditions
  // -----------------------------
  void getInitConditions( TacsScalar v[], TacsScalar dv[],
                          TacsScalar ddv[], const TacsScalar X[] ) override;

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time,
                    TacsScalar res[],
                    const TacsScalar X[], const TacsScalar v[],
                    const TacsScalar dv[], const TacsScalar ddv[] ) override;

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
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
  
  void addAdjResProduct( double time,
                         double scale,
                         TacsScalar dfdx[],
                         int dvLen,
                         const TacsScalar psi[],
                         const TacsScalar Xpts[],
                         const TacsScalar v[],
                         const TacsScalar dv[],
                         const TacsScalar ddv[] ) override;
 
  void addAdjResXptProduct( double time, double scale,
                            TacsScalar dfdx[],
                            const TacsScalar psi[],
                            const TacsScalar Xpts[],
                            const TacsScalar v[],
                            const TacsScalar dv[],
                            const TacsScalar ddv[] ) override;
 
  void setDesignVars( const TacsScalar dvs[], int numDVs ) override{
    delem->setDesignVars(dvs, numDVs);
  }
  void getDesignVars( TacsScalar dvs[], int numDVs ) override{
    delem->getDesignVars(dvs, numDVs);
  }
  void getDesignVarRange( TacsScalar lowerBound[],
                          TacsScalar upperBound[],
                          int numDVs ) override{
    delem->getDesignVarRange(lowerBound, upperBound, numDVs);
  }
 
  // Invoke this function to update this element through user supplied callback
  //---------------------------------------------------------------------------
  void updateElement(TACSElement* elem, TacsScalar* vals){
    if (this->update && pyptr){
      this->update(elem, vals, pyptr);
    } else {
      if (this->update) this->update(elem, vals, NULL);
    }
  }

  TACSElement* getDeterministicElement(){
    return this->delem;
  };
  
  // Callback function to update the parameters of element
  void (*update)(TACSElement*, TacsScalar*, void*);
  PyObject *pyptr; 

 protected:
  TACSElement *delem;
  ParameterContainer *pc;

 private:
  // Stochastic element information
  int num_nodes;
  int vars_per_node;
};

#endif
