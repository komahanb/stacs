#include "TACSStochasticElement.h"

namespace{
  /*
    Predetermine whether the jacobian block is nonzero
  */
  bool nonzero(const int nvars, int *dmapi, int *dmapj, int *dmapk){
    int filter[nvars];
    for( int i = 0; i < nvars; i++ ){
      if (abs(dmapi[i] - dmapj[i]) <= dmapk[i]){
        filter[i] = 1;
      } else {
        filter[i] = 0;
      }
    }
    int prod = 1;
    for( int i = 0; i < nvars; i++ ){
      prod *= filter[i];
    }
    return prod;
  }

  /*
    Place entry into the matrix location
  */
  TacsScalar getElement(TacsScalar *A, int size, int row, int col) {
    return A[size * row + col];
  };

  /*
    Add entry into the matrix location
  */
  void addElement(TacsScalar *J, int size, int row, int col, TacsScalar value) {
    J[size * row + col] += value;
  };

  /*
    Display the matrix sparsity pattern
  */
  void printSparsity( TacsScalar *J, int size,  double tol = 1.0e-12 ){
    for (int ii = 0; ii < size; ii++){
      printf("%2d ", ii);
      for (int jj = 0; jj < size; jj++){
        if (abs(getElement(J, size, ii, jj)) > tol) {
          printf("%2s", "x");
        } else {
          printf("%2s", " ");
        }
      }
      printf("\n");
    }
  }

  void getDeterministicAdjoint( ParameterContainer *pc, 
                                TACSElement *delem,
                                TACSElement *selem, 
                                const TacsScalar v[],
                                TacsScalar *zq,
                                TacsScalar *uq
                                ){
    int ndvpn   = delem->getVarsPerNode();
    int nsvpn   = selem->getVarsPerNode();
    int nddof   = delem->getNumVariables();
    int nsdof   = selem->getNumVariables();
    int nsterms = pc->getNumBasisTerms();
    int nnodes  = selem->getNumNodes();
    
    memset(uq  , 0, nddof*sizeof(TacsScalar));

    // Evaluate the basis at quadrature node and form the state
    // vectors
    for (int n = 0; n < nnodes; n++){
      for (int k = 0; k < nsterms; k++){
        TacsScalar psikz = pc->basis(k,zq);
        int lptr = n*ndvpn;
        int gptr = n*nsvpn + k*ndvpn;
        for (int d = 0; d < ndvpn; d++){        
          uq[lptr+d] += v[gptr+d]*psikz;
        }
      }
    }
  } 

  void getDeterministicStates( ParameterContainer *pc, 
                               TACSElement *delem,
                               TACSElement *selem, 
                               const TacsScalar v[],
                               const TacsScalar dv[],
                               const TacsScalar ddv[], 
                               TacsScalar *zq,
                               TacsScalar *uq,
                               TacsScalar *udq,
                               TacsScalar *uddq
                               ){
    int ndvpn   = delem->getVarsPerNode();
    int nsvpn   = selem->getVarsPerNode();
    int nddof   = delem->getNumVariables();
    int nsdof   = selem->getNumVariables();
    int nsterms = pc->getNumBasisTerms();
    int nnodes  = selem->getNumNodes();

    memset(uq  , 0, nddof*sizeof(TacsScalar));
    memset(udq , 0, nddof*sizeof(TacsScalar));
    memset(uddq, 0, nddof*sizeof(TacsScalar));

    // Evaluate the basis at quadrature node and form the state
    // vectors
    for (int n = 0; n < nnodes; n++){
      for (int k = 0; k < nsterms; k++){
        TacsScalar psikz = pc->basis(k,zq);
        int lptr = n*ndvpn;
        int gptr = n*nsvpn + k*ndvpn;
        for (int d = 0; d < ndvpn; d++){        
          uq[lptr+d] += v[gptr+d]*psikz;
          udq[lptr+d] += dv[gptr+d]*psikz;
          uddq[lptr+d] += ddv[gptr+d]*psikz;
        }
      }
    }
  } 

 
}

TACSStochasticElement::TACSStochasticElement( TACSElement *_delem,
                                              ParameterContainer *_pc,
                                              void (*_update)(TACSElement*, TacsScalar*, void*) ){
  // Store the deterministic element
  delem = _delem;
  delem->incref();

  // Set callback for element update
  update = _update;

  // Set the component numner of this element
  setComponentNum(delem->getComponentNum());

  // Parameter container
  pc = _pc;

  // Set number of dofs
  num_nodes     = delem->getNumNodes();
  vars_per_node = pc->getNumBasisTerms()*delem->getVarsPerNode();
}

TACSStochasticElement::~TACSStochasticElement(){
  this->delem->decref();
  this->delem = NULL;
  this->pc = NULL;
}

/*
  TACS Element member functions
*/
int TACSStochasticElement::getVarsPerNode() {
  return vars_per_node;
}

int TACSStochasticElement::getNumNodes() {
  return num_nodes;
}

/*
  Return the Initial conditions after projection
*/
void TACSStochasticElement::getInitConditions( int elemIndex,
                                               const TacsScalar X[],
                                               TacsScalar v[],
                                               TacsScalar dv[],
                                               TacsScalar ddv[] ){
  const int ndvpn   = delem->getVarsPerNode();
  const int nsvpn   = this->getVarsPerNode();
  const int nddof   = delem->getNumVariables();
  const int nsdof   = this->getNumVariables();
  const int nsterms = pc->getNumBasisTerms();
  const int nnodes  = this->getNumNodes();

  // Space for quadrature points and weights
  const int nsparams = pc->getNumParameters();
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;

  // Create space for states
  TacsScalar *uq    = new TacsScalar[nddof];
  TacsScalar *udq   = new TacsScalar[nddof];
  TacsScalar *uddq  = new TacsScalar[nddof];

  TacsScalar *utmpk    = new TacsScalar[nddof];
  TacsScalar *udtmpk   = new TacsScalar[nddof];
  TacsScalar *uddtmpk  = new TacsScalar[nddof];

  const int nqpts = pc->getNumQuadraturePoints();

  //  Projection of initial conditions and return
  for (int k = 0; k < nsterms; k++){

    memset(utmpk  , 0, nddof*sizeof(TacsScalar));
    memset(udtmpk , 0, nddof*sizeof(TacsScalar));
    memset(uddtmpk, 0, nddof*sizeof(TacsScalar));

    for (int q = 0; q < nqpts; q++){
      // Get the quadrature points and weights
      wq = pc->quadrature(q, zq, yq);

      // Set the parameter values into the element
      updateElement(delem, yq);

      // reset the states to zero
      memset(uq  , 0, nddof*sizeof(TacsScalar));
      memset(udq , 0, nddof*sizeof(TacsScalar));
      memset(uddq, 0, nddof*sizeof(TacsScalar));

      // Fetch the deterministic element residual
      delem->getInitConditions(elemIndex, X, uq, udq, uddq);

      // Project the determinic states onto the stochastic basis and
      // place in global state array
      TacsScalar scale = pc->basis(k,zq)*wq;
      for (int c = 0; c < nddof; c++){
        utmpk[c] += uq[c]*scale;
        udtmpk[c] += udq[c]*scale;
        uddtmpk[c] += uddq[c]*scale;
      }

    } // quadrature

    // Store the initial conditions in termwise order
    for (int n = 0; n < nnodes; n++){
      int lptr = n*ndvpn;
      int gptr = n*nsvpn + k*ndvpn;
      for (int d = 0; d < ndvpn; d++){        
        v[gptr+d] = utmpk[lptr+d];
        dv[gptr+d] = udtmpk[lptr+d];
        ddv[gptr+d] = uddtmpk[lptr+d];
      }
    }
  }

  // clear the heap
  delete [] uq;
  delete [] udq;
  delete [] uddq;
  delete [] utmpk;
  delete [] udtmpk;
  delete [] uddtmpk;
  delete [] zq;
  delete [] yq;
}

/*
  Compute the residual of the governing equations
*/
void TACSStochasticElement::addResidual( int elemIndex,
                                         double time,
                                         const TacsScalar X[],
                                         const TacsScalar v[],
                                         const TacsScalar dv[],
                                         const TacsScalar ddv[],
                                         TacsScalar res[] ){
  const int ndvpn   = delem->getVarsPerNode();
  const int nsvpn   = this->getVarsPerNode();
  const int nddof   = delem->getNumVariables();
  const int nsdof   = this->getNumVariables();
  const int nsterms = pc->getNumBasisTerms();
  const int nnodes  = this->getNumNodes();

  // Space for quadrature points and weights
  const int nsparams = pc->getNumParameters();
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;

  // Create space for fetching deterministic residuals and states
  TacsScalar *uq    = new TacsScalar[nddof];
  TacsScalar *udq   = new TacsScalar[nddof];
  TacsScalar *uddq  = new TacsScalar[nddof];
  TacsScalar *resq  = new TacsScalar[nddof];
  TacsScalar *rtmpi = new TacsScalar[nddof];

  const int nqpts = pc->getNumQuadraturePoints();

  for (int i = 0; i < nsterms; i++){

    memset(rtmpi, 0, nddof*sizeof(TacsScalar));

    for (int q = 0; q < nqpts; q++){

      // Get the quadrature points and weights
      wq = pc->quadrature(q, zq, yq);

      // Set the parameter values into the element
      updateElement(delem, yq);

      // reset the states and residuals
      memset(resq, 0, nddof*sizeof(TacsScalar));

      // Evaluate the basis at quadrature node and form the state
      // vectors
      getDeterministicStates(pc, delem, this, v, dv, ddv, zq, 
                             uq, udq, uddq);

      // Fetch the deterministic element residual
      delem->addResidual(elemIndex, time, X, uq, udq, uddq, resq);

      //  Project the determinic element residual onto the
      //  stochastic basis and place in global residual array
      TacsScalar scale = pc->basis(i,zq)*wq;
      for (int c = 0; c < nddof; c++){
        rtmpi[c] += resq[c]*scale;
      }

    } // quadrature

    // Store i-th projected Residual into stochastic array
    for (int n = 0; n < nnodes; n++){
      int lptr = n*ndvpn;
      int gptr = n*nsvpn + i*ndvpn;
      for (int d = 0; d < ndvpn; d++){        
        res[gptr+d] += rtmpi[lptr+d];
      }
    }

  } // end nsterms

  // clear the heap
  delete [] rtmpi;
  delete [] resq;
  delete [] uq;
  delete [] udq;
  delete [] uddq;
  delete [] zq;
  delete [] yq;
}

void TACSStochasticElement::addJacobian( int elemIndex,
                                         double time,
                                         TacsScalar alpha,
                                         TacsScalar beta,
                                         TacsScalar gamma,
                                         const TacsScalar X[],
                                         const TacsScalar v[],
                                         const TacsScalar dv[],
                                         const TacsScalar ddv[],
                                         TacsScalar res[],
                                         TacsScalar mat[] ){
  // Call the residual implementation
  addResidual(elemIndex, time, X, v, dv, ddv, res);

  const int ndvpn   = delem->getVarsPerNode();
  const int nsvpn   = this->getVarsPerNode();
  const int nddof   = delem->getNumVariables();
  const int nsdof   = this->getNumVariables();
  const int nsterms = pc->getNumBasisTerms();
  const int nnodes  = this->getNumNodes();

  // Space for quadrature points and weights
  const int nsparams = pc->getNumParameters();
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;

  // polynomial degrees
  int dmapi[nsparams], dmapj[nsparams], dmapf[nsparams];

  // Assume all parameters to be of degree atmost 3 (faster)
  for (int i = 0; i <nsparams; i++){
    dmapf[i] = 3;
  }

  // Create space for fetching deterministic residuals and states
  TacsScalar *uq    = new TacsScalar[nddof];
  TacsScalar *udq   = new TacsScalar[nddof];
  TacsScalar *uddq  = new TacsScalar[nddof];
  TacsScalar *A     = new TacsScalar[nddof*nddof];
  TacsScalar *resq  = new TacsScalar[nddof];

  const int nqpts = pc->getNumQuadraturePoints();

  for (int i = 0; i < nsterms; i++){

    pc->getBasisParamDeg(i, dmapi);

    for (int j = 0; j < nsterms; j++){

      pc->getBasisParamDeg(j, dmapj);

      if (1){ // nonzero(nsparams, dmapi, dmapj, dmapf)){

        memset(A, 0, nddof*nddof*sizeof(TacsScalar));

        for (int q = 0; q < nqpts; q++){

          // Get quadrature points
          wq = pc->quadrature(q, zq, yq);

          // Set the parameter values into the element
          this->updateElement(this->delem, yq);

          // Evaluate the basis at quadrature node and form the state
          // vectors
          getDeterministicStates(pc, delem, this, v, dv, ddv, zq, 
                                 uq, udq, uddq);


          // Fetch the deterministic element residual
          TacsScalar scale = pc->basis(i,zq)*pc->basis(j,zq)*wq;
          this->delem->addJacobian(elemIndex,
                                   time,
                                   scale*alpha, scale*beta, scale*gamma,
                                   X, uq, udq, uddq,
                                   resq,
                                   A);
        } // quadrature

        // Place the (i,j)-projected block into the stochastic block
        for (int ni = 0; ni < nnodes; ni++){
          int liptr = ni*ndvpn;
          int giptr = ni*nsvpn + i*ndvpn;
          for (int di = 0; di < ndvpn; di++){
            for (int nj = 0; nj < nnodes; nj++){
              int ljptr = nj*ndvpn;
              int gjptr = nj*nsvpn + j*ndvpn;
              for (int dj = 0; dj < ndvpn; dj++){
                addElement(mat, nsdof, 
                           giptr + di, gjptr + dj,
                           getElement(A, nddof, 
                                      liptr + di, ljptr + dj));
              }
            }
          }
        }

      } // nonzero

    } // end j

  } // end i

  //  printSparsity(mat, nddof*nsterms);

  // clear the heap
  delete [] A;
  delete [] resq;
  delete [] uq;
  delete [] udq;
  delete [] uddq;
  delete [] zq;
  delete [] yq;
}

int TACSStochasticElement::evalPointQuantity( int elemIndex, int quantityType, double time,
                                              int N, double pt[], const TacsScalar Xpts[],
                                              const TacsScalar v[], const TacsScalar dv[],
                                              const TacsScalar ddv[], TacsScalar *quantity ) {
  const int ndvpn   = delem->getVarsPerNode();
  const int nsvpn   = this->getVarsPerNode();
  const int nddof   = delem->getNumVariables();
  const int nnodes  = this->getNumNodes();
  const int nsterms = pc->getNumBasisTerms();

  // Space for quadrature points and weights
  const int nsparams = pc->getNumParameters();
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;
  
  // Create space for deterministic states at each quadrature node in y
  TacsScalar *uq     = new TacsScalar[nddof];
  TacsScalar *udq    = new TacsScalar[nddof];
  TacsScalar *uddq   = new TacsScalar[nddof];
  
  // Space to project each function in stochastic space and store
  // const  int ndquants = this->delem->getNumPointQuantities();
  const int ndquants = this->delem->evalPointQuantity(elemIndex,
                                                      quantityType,
                                                      time, N, pt, Xpts, 
                                                      v, dv, ddv, 
                                                      uq); // dummy 
  const int nsquants = nsterms*ndquants;

  TacsScalar *ftmpq  = new TacsScalar[ndquants];
  TacsScalar *ftmpi  = new TacsScalar[ndquants];  

  const int nqpts = pc->getNumQuadraturePoints();
  
  for (int i = 0; i < nsterms; i++){

    memset(ftmpi, 0, ndquants*sizeof(TacsScalar));

    for (int q = 0; q < nqpts; q++){

      // Get the quadrature points and weights
      wq = pc->quadrature(q, zq, yq);

      // Set the parameter values into the element
      this->updateElement(delem, yq);

      // reset the states and residuals
      memset(ftmpq, 0, ndquants*sizeof(TacsScalar));

      // Evaluate the basis at quadrature node and form the state
      // vectors
      getDeterministicStates(pc, delem, this, v, dv, ddv, zq, 
                             uq, udq, uddq);

      // Fetch the deterministic element residual
      int count = this->delem->evalPointQuantity(elemIndex,
                                                 quantityType,
                                                 time, N, pt,
                                                 Xpts, uq, udq, uddq,
                                                 ftmpq);
      
      // Project the determinic quantities onto the stochastic basis
      // and place in stochastic function array
      TacsScalar scale = pc->basis(i,zq)*wq;
      for (int c = 0; c < ndquants; c++){
        ftmpi[c] += ftmpq[c]*scale;
      }

    } // quadrature

    // Store i-th projected Residual into stochastic array
    for (int d = 0; d < ndquants; d++){        
      quantity[d*nsterms+i] = ftmpi[d];
      // printf("nsterms = %d ndquants = %d, i = %d local [%d]  global[%d] \n", nsterms, ndquants,
      // i, d, d*nsterms + i);
    }

  } // nterms

  delete [] zq;
  delete [] yq;

  delete [] uq;
  delete [] udq;
  delete [] uddq;

  delete [] ftmpq;
  delete [] ftmpi;

  return nsquants;
}

/*
  Stochastic Adjoint residual product implementation
*/
void TACSStochasticElement::addAdjResProduct( int elemIndex, double time,
                                              TacsScalar scale,
                                              const TacsScalar psi[],
                                              const TacsScalar Xpts[],
                                              const TacsScalar v[],
                                              const TacsScalar dv[],
                                              const TacsScalar ddv[],
                                              int dvLen,
                                              TacsScalar dfdx[] ){

  //  printf("TACSStochasticElement::addAdjResProduct \n");

  const int ndvpn   = delem->getVarsPerNode();
  const int nsvpn   = this->getVarsPerNode();
  const int nddof   = delem->getNumVariables();
  const int nnodes  = this->getNumNodes();
  const int nsterms = pc->getNumBasisTerms();

  // Space for quadrature points and weights
  const int nsparams = pc->getNumParameters();
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;
  
  // Create space for deterministic states at each quadrature node in y
  TacsScalar *uq     = new TacsScalar[nddof];
  TacsScalar *udq    = new TacsScalar[nddof];
  TacsScalar *uddq   = new TacsScalar[nddof];
  TacsScalar *psiq   = new TacsScalar[nddof];
  TacsScalar *dfdxj  = new TacsScalar[dvLen]; // check if this is one function at a time

  const int nqpts = pc->getNumQuadraturePoints();
  
  for (int j = 0; j < 1; j++){

    memset(dfdxj, 0, dvLen*sizeof(TacsScalar));

    for (int q = 0; q < nqpts; q++){

      // Get the quadrature points and weights
      wq = pc->quadrature(q, zq, yq);

      TacsScalar wt = pc->basis(j,zq)*wq;

      // Set the parameter values into the element
      this->updateElement(delem, yq);

      // Form deterministi states and adjoint vectors from global array
      getDeterministicStates(pc, delem, this, v, dv, ddv, zq, uq, udq, uddq);
      getDeterministicAdjoint(pc, delem, this, psi, zq, psiq);

      delem->addAdjResProduct(elemIndex, time, wt*scale,
                              psiq, Xpts, uq, udq, uddq,
                              dvLen, dfdxj);

    } // end quadrature

    // need to be careful with nodewise placement of dvs
    for (int n = 0; n < dvLen; n++){
      dfdx[n] += dfdxj[n];
    }
    
    // should we store design variablewise or stochastic termwise
    /*
    // Termwise storage of dfdxj in to dfdx
    for (int i = 0; i < dvLen; i++){
    dfdx[j*nsterms+i] +=  dfdxj[i];
    }
    
    // fix
    // Design variablewise storage of projected dfdx
    for (int i = 0; i < dvLen; i++){
      dfdx[i*nsterms+j] +=  dfdxj[i];
    }
    */
  } // nsterms

  // fix
  // how is 'dfdx' bigger for stochastic element? are we multiplying
  // the nsterms with num det design variables?

  // free up heap
  delete [] zq;
  delete [] yq;
  delete [] uq;
  delete [] udq;
  delete [] uddq;
  delete [] psiq;
  delete [] dfdxj;
}
