#include "TACSAssembler.h"
#include "TACSKSStochasticFMeanFunction.h"
#include "TACSStochasticElement.h"

namespace {

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

TACSKSStochasticFMeanFunction::TACSKSStochasticFMeanFunction( TACSAssembler *tacs,
                                                              TACSFunction *dfunc,
                                                              ParameterContainer *pc,
                                                              int quantityType,
                                                              int moment_type,
                                                              double ksWeight )
  : TACSFunction(dfunc->getAssembler(), 
                 dfunc->getDomainType(), 
                 dfunc->getStageType(),
                 0)
{ 
  this->tacs_comm = tacs->getMPIComm();
  this->dfunc = dfunc;
  this->dfunc->incref();
  this->pc = pc;
  this->quantityType = quantityType;
  this->moment_type = moment_type;
  this->ksWeight = ksWeight;
  this->nsqpts  = pc->getNumQuadraturePoints();
  this->nsterms = pc->getNumBasisTerms();
  this->fvals    = new TacsScalar[nsterms*nsqpts];
  this->ksSum    = new TacsScalar[nsterms*nsqpts];
  this->maxValue = new TacsScalar[nsterms*nsqpts];
}

TACSKSStochasticFMeanFunction::~TACSKSStochasticFMeanFunction()
{
  delete [] this->ksSum;
  delete [] this->maxValue;
  delete [] this->fvals;
  delete [] this->ksSum;
  delete [] this->maxValue;
}

void TACSKSStochasticFMeanFunction::initEvaluation( EvaluationType ftype )
{
  if (ftype == TACSFunction::INITIALIZE){
    for (int k = 0; k < nsterms*nsqpts; k++){
      this->maxValue[k] = -1.0e20;
    }     
  }
  else if (ftype == TACSFunction::INTEGRATE){
    for (int k = 0; k < nsterms*nsqpts; k++){
      this->ksSum[k] = 0.0;
    }     
  }
}

void TACSKSStochasticFMeanFunction::elementWiseEval( EvaluationType evalType,
                                                     int elemIndex,
                                                     TACSElement *element,
                                                     double time,
                                                     TacsScalar tscale,
                                                     const TacsScalar Xpts[],
                                                     const TacsScalar v[],
                                                     const TacsScalar dv[],
                                                     const TacsScalar ddv[] )
{
  //printf("TACSStochasticVarianceFunction::elementWiseEval %d\n", elemIndex);
  TACSStochasticElement *selem = dynamic_cast<TACSStochasticElement*>(element);
  if (!selem) {
    printf("Casting to stochastic element failed; skipping elemenwiseEval");
  };
  
  TACSElement *delem = selem->getDeterministicElement();
  const int nqpts    = pc->getNumQuadraturePoints();
  const int nsparams = pc->getNumParameters();
  const int ndvpn    = delem->getVarsPerNode();
  const int nsvpn    = selem->getVarsPerNode();
  const int nddof    = delem->getNumVariables();
  const int nnodes   = selem->getNumNodes();  
  
  // Space for quadrature points and weights
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;
  
  // Create space for deterministic states at each quadrature node in y
  TacsScalar *uq     = new TacsScalar[nddof];
  TacsScalar *udq    = new TacsScalar[nddof];
  TacsScalar *uddq   = new TacsScalar[nddof];

  for (int j = 0; j < nsterms; j++){

    // Stochastic Integration
    for (int q = 0; q < nqpts; q++){

      // Get the quadrature points and weights for mean
      wq = pc->quadrature(q, zq, yq);

      // Set the parameter values into the element
      selem->updateElement(delem, yq);

      // Form the state vectors
      getDeterministicStates(pc, delem, selem, 
                             v, dv, ddv, zq, 
                             uq, udq, uddq);

      {
        // Get the number of quadrature points for this delem
        const int numGauss = 1; //delem->getNumGaussPts();
        const int numDisps = delem->getNumVariables();
        const int numNodes = delem->getNumNodes();

        for ( int i = 0; i < numGauss; i++ ){
      
          // Get the Gauss points one at a time
          TacsScalar weight = 1.0; //delem->getGaussWtsPts(i, pt);
          double pt[3] = {0.0,0.0,0.0};
          const int N = 1;
          //  delem->getShapeFunctions(pt, ctx->N);
   
          // Evaluate the dot-product with the displacements
          // const double *N = ctx->N;

          TacsScalar value = 0.0;
          delem->evalPointQuantity(elemIndex,
                                   this->quantityType,
                                   time, N, pt,
                                   Xpts, uq, udq, uddq,
                                   &value);

          if (evalType == TACSFunction::INITIALIZE){      
            // Reset maxvalue if needed
            if (TacsRealPart(value) > TacsRealPart(maxValue[j*nsqpts+q])){
              maxValue[j*nsqpts+q] = value;
            }      
          } else {
            // Add up the contribution from the quadrature
            // delem->getDetJacobian(pt, Xpts);
            TacsScalar h = 1.0;
            ksSum[j*nsqpts+q] += tscale*exp(ksWeight*(value - maxValue[j*nsqpts+q]));
          }      

        } // spatial integration

      }
      
    } // end yloop

  } // nsterms

  // clear allocated heap
  delete [] zq;
  delete [] yq;
  delete [] uq;
  delete [] udq;
  delete [] uddq;
}

void TACSKSStochasticFMeanFunction::finalEvaluation( EvaluationType evalType )
{
  if (evalType == TACSFunction::INITIALIZE){
    TacsScalar temp;
    for (int q = 0; q < nsterms*nsqpts; q++){
      temp = maxValue[q];
      MPI_Allreduce(&temp, &maxValue[q], 1, TACS_MPI_TYPE, TACS_MPI_MAX, this->tacs_comm);
    }
  } else {
    TacsScalar temp;
    for (int q = 0; q < nsterms*nsqpts; q++){
      temp = ksSum[q];
      MPI_Allreduce(&temp, &ksSum[q], 1, TACS_MPI_TYPE, MPI_SUM, this->tacs_comm);
    }
    for (int k = 0; k < nsterms; k++){
      for (int q = 0; q < nsqpts; q++){
        fvals[k*nsqpts+q] = maxValue[k*nsqpts+q] + log(ksSum[k*nsqpts+q])/ksWeight;
      }
    }

  }
}

/**
   Get the value of the function
*/
TacsScalar TACSKSStochasticFMeanFunction::getFunctionValue(){  
  return getExpectation();
}

TacsScalar TACSKSStochasticFMeanFunction::getExpectation(){
  // Finish up stochastic integration
  const int nsparams = pc->getNumParameters();
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;
  TacsScalar fmean = 0.0;    
  for (int q = 0; q < nsqpts; q++){
    TacsScalar wq = pc->quadrature(q, zq, yq);
    fmean += wq*pc->basis(0,zq)*fvals[0*nsqpts+q];
  }
  delete [] zq;
  delete [] yq;
  return fmean;
}
 
TacsScalar TACSKSStochasticFMeanFunction::getVariance(){
  TacsScalar fvar = 0.0;
  for (int k = 1; k < nsterms; k++){
    fvar += fvals[k]*fvals[k];
  }
  return fvar;
}

void TACSKSStochasticFMeanFunction::getElementSVSens( int elemIndex, TACSElement *element,
                                                      double time,
                                                      TacsScalar alpha, TacsScalar beta,
                                                      TacsScalar gamma,
                                                      const TacsScalar Xpts[],
                                                      const TacsScalar v[],
                                                      const TacsScalar dv[],
                                                      const TacsScalar ddv[],
                                                      TacsScalar dfdu[] ){
  if (RealPart(ksSum[0]) < 1.0e-15){
    printf("Error: Evaluate the functions before derivatives \n");
  }

  TACSStochasticElement *selem = dynamic_cast<TACSStochasticElement*>(element);
  if (!selem) {
    printf("Casting to stochastic element failed; skipping elemenwiseEval");
  };

  
  TACSElement *delem = selem->getDeterministicElement();
  const int nsterms  = pc->getNumBasisTerms();
  const int nqpts    = pc->getNumQuadraturePoints();
  const int nsparams = pc->getNumParameters();
  const int ndvpn    = delem->getVarsPerNode();
  const int nsvpn    = selem->getVarsPerNode();
  const int nddof    = delem->getNumVariables();
  const int nnodes   = selem->getNumNodes();  
  const int nsdof    = selem->getNumVariables();
  memset(dfdu, 0, nsdof*sizeof(TacsScalar));

  // j-th project
  TacsScalar *dfduj  = new TacsScalar[nddof];  
  
  // Space for quadrature points and weights
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;
  
  // Create space for deterministic states at each quadrature node in y
  TacsScalar *uq     = new TacsScalar[nddof];
  TacsScalar *udq    = new TacsScalar[nddof];
  TacsScalar *uddq   = new TacsScalar[nddof];

  for (int j = 0; j < nsterms; j++){

    memset(dfduj, 0, nddof*sizeof(TacsScalar));
    
    // Stochastic Integration
    for (int q = 0; q < nsqpts; q++){

      // Get the quadrature points and weights for mean
      wq = pc->quadrature(q, zq, yq);
      TacsScalar wt = pc->basis(j,zq)*wq;
    
      // Set the parameter values into the element
      selem->updateElement(delem, yq);

      // Form the state vectors
      getDeterministicStates(pc, delem, selem, 
                             v, dv, ddv, zq, 
                             uq, udq, uddq);

      { 

        // Get the number of quadrature points for this element
        const int numGauss = 1; //delem->getNumGaussPts();
        const int numDisps = delem->getNumVariables();
        const int numNodes = delem->getNumNodes();
    
        for ( int i = 0; i < numGauss; i++ ){      
          TacsScalar weight       = 1.0; //delem->getGaussWtsPts(i, pt);
          double pt[3]        = {0.0,0.0,0.0};
          const int N         = 1;

          TacsScalar quantity = 0.0;
          delem->evalPointQuantity(elemIndex,
                                   this->quantityType,
                                   time, N, pt,
                                   Xpts, uq, udq, uddq,
                                   &quantity);        
          
          TacsScalar ksPtWeight = 0.0;
          ksPtWeight = exp(ksWeight*(quantity - maxValue[j*nsqpts+q]))/ksSum[j*nsqpts+q];
          // ksPtWeight *= weight*detJ;

          //          printf("%.17e %.17e quantity = %.17e ksptweight = %.17e \n", maxValue[j*nsqpts+q], ksSum[j*nsqpts+q], quantity, ksPtWeight);
          TacsScalar dfdq = ksPtWeight;
          delem->addPointQuantitySVSens(elemIndex,
                                        this->quantityType,
                                        time,
                                        wt*alpha*ksPtWeight,
                                        wt*beta*ksPtWeight,
                                        wt*gamma*ksPtWeight,
                                        N, pt,
                                        Xpts, uq, udq, uddq,
                                        &dfdq, dfduj);
        } // spatial integration

      }

    } // probabilistic integration
    
    // Store j-th projected sv sens into stochastic array
    for (int n = 0; n < nnodes; n++){
      int lptr = n*ndvpn;
      int gptr = n*nsvpn + j*ndvpn;
      for (int d = 0; d < ndvpn; d++){        
        dfdu[gptr+d] = dfduj[lptr+d];
      }
    }

  } // end nsterms

  // clear allocated heap
  delete [] dfduj;
  delete [] zq;
  delete [] yq;
  delete [] uq;
  delete [] udq;
  delete [] uddq; 
}

void TACSKSStochasticFMeanFunction::addElementDVSens( int elemIndex, TACSElement *element,
                                                      double time, TacsScalar scale,
                                                      const TacsScalar Xpts[],
                                                      const TacsScalar v[],
                                                      const TacsScalar dv[],
                                                      const TacsScalar ddv[],
                                                      int dvLen,
                                                      TacsScalar dfdx[] ){


  TACSStochasticElement *selem = dynamic_cast<TACSStochasticElement*>(element);
  if (!selem) {
    printf("Casting to stochastic element failed; skipping elemenwiseEval");
  };
  TACSElement *delem  = selem->getDeterministicElement();
  
  const int nsterms   = pc->getNumBasisTerms();
  const int nsparams  = pc->getNumParameters();
  const int ndvpn     = delem->getVarsPerNode();
  const int nsvpn     = selem->getVarsPerNode();
  const int nddof     = delem->getNumVariables();
  const int nnodes    = selem->getNumNodes();  
  const int dvpernode = delem->getDesignVarsPerNode();

  // j-th projection of dfdx array
  TacsScalar *dfdxj  = new TacsScalar[dvLen];
  
  // Space for quadrature points and weights
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;
  
  // Create space for deterministic states at each quadrature node in y
  TacsScalar *uq     = new TacsScalar[nddof];
  TacsScalar *udq    = new TacsScalar[nddof];
  TacsScalar *uddq   = new TacsScalar[nddof];
  
  // int nterms;
  // if (moment_type == 0){
  //   nterms = 1;
  // } else {
  //   nterms = 
  // }

  for (int j = 0; j < 1; j++){   // not sure, other contributions are zero, so this does not affect

    memset(dfdxj, 0, dvLen*sizeof(TacsScalar));
    
    // Stochastic Integration
    for (int q = 0; q < nsqpts; q++){

      // Get the quadrature points and weights for mean
      wq = pc->quadrature(q, zq, yq);
      TacsScalar wt = pc->basis(j,zq)*wq;
    
      // Set the parameter values into the element
      selem->updateElement(delem, yq);

      // form deterministic states      
      getDeterministicStates(pc, delem, selem, v, dv, ddv, zq, uq, udq, uddq);


      { 

        // Get the number of quadrature points for this element
        const int numGauss = 1; //delem->getNumGaussPts();
        const int numDisps = delem->getNumVariables();
        const int numNodes = delem->getNumNodes();
    
        for ( int i = 0; i < numGauss; i++ ){      
          TacsScalar weight       = 1.0; //delem->getGaussWtsPts(i, pt);
          double pt[3]        = {0.0,0.0,0.0};
          const int N         = 1;

          TacsScalar quantity = 0.0;
          delem->evalPointQuantity(elemIndex,
                                   this->quantityType,
                                   time, N, pt,
                                   Xpts, uq, udq, uddq,
                                   &quantity);        
          
          TacsScalar ksPtWeight = 0.0;
          ksPtWeight = exp(ksWeight*(quantity - maxValue[j*nsqpts+q]))/ksSum[j*nsqpts+q];
          // ksPtWeight *= weight*detJ;

          //printf("%.17e %.17e ksptweight = %.17e \n", maxValue[j*nsqpts+q], ksSum[j*nsqpts+q], ksPtWeight);

          // Call the underlying element and get the design variable sensitivities
          TacsScalar _dfdq = ksPtWeight; 
          delem->addPointQuantityDVSens( elemIndex, 
                                         this->quantityType,
                                         time, ksPtWeight*wt*scale,
                                         N, pt,
                                         Xpts, uq, udq, uddq, &_dfdq, 
                                         dvLen, dfdxj ); 

        } // spatial integration
        
      }
      
    } // end yloop

    // for (int n = 0; n < nnodes; n++){
    //   int lptr = n*dvpernode;
    //   int gptr = n*dvpernode + j*dvpernode;
    //   for (int i = 0; i < dvpernode; i++){
    //     dfdx[gptr+n] += dfdxj[lptr+n];
    //   }            
    // }

    // need to be careful with nodewise placement of dvsx
    for (int n = 0; n < dvLen; n++){
      //printf("term %d dfdx[%d] = %.17e %.17e \n", j, n, dfdx[n], dfdxj[n]);
      dfdx[n] += dfdxj[n];
    }
    
  } // end nsterms

  // clear allocated heap
  delete [] zq;
  delete [] yq;
  delete [] uq;
  delete [] udq;
  delete [] uddq;
  delete [] dfdxj;
}

