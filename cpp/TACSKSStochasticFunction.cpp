#include "TACSAssembler.h"
#include "TACSKSStochasticFunction.h"
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

TACSKSStochasticFunction::TACSKSStochasticFunction( TACSAssembler *tacs,
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

TACSKSStochasticFunction::~TACSKSStochasticFunction()
{
  delete [] this->fvals;
  delete [] this->ksSum;
  delete [] this->maxValue;
}

void TACSKSStochasticFunction::initEvaluation( EvaluationType ftype )
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

void TACSKSStochasticFunction::elementWiseEval( EvaluationType evalType,
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
        TACSElementBasis *basis = delem->getElementBasis();

        if (basis){

          for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){

            double pt[3];
            double weight = basis->getQuadraturePoint(i, pt);
            TacsScalar value = 0.0;
            int count = delem->evalPointQuantity(elemIndex,
                                                 this->quantityType,
                                                 time, i, pt,
                                                 Xpts, uq, udq, uddq,
                                                 &value);
          
            if (evalType == TACSFunction::INITIALIZE){      
              // Reset maxvalue if needed
              if (TacsRealPart(value) > TacsRealPart(maxValue[j*nsqpts+q])){
                maxValue[j*nsqpts+q] = value;
              }      
            } else {
              // Evaluate the determinant of the Jacobian
              TacsScalar Xd[9], J[9];
              TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);
              ksSum[j*nsqpts+q] += tscale*weight*detJ*exp(ksWeight*(value - maxValue[j*nsqpts+q]));
            }      

          } // spatial integration

        } // basis

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

void TACSKSStochasticFunction::finalEvaluation( EvaluationType evalType )
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
TacsScalar TACSKSStochasticFunction::getFunctionValue(){ 
  // Finish up stochastic integration
  const int nsparams = pc->getNumParameters();
  TacsScalar *zq = new TacsScalar[nsparams];
  TacsScalar *yq = new TacsScalar[nsparams];
  TacsScalar wq;
  
  TacsScalar fmean = 0.0;
  for (int q = 0; q < nsqpts; q++){
    wq = pc->quadrature(q, zq, yq);
    fmean += wq*pc->basis(0,zq)*fvals[0*nsqpts+q];
  }

  TacsScalar ffmean = 0.0;       
  if (moment_type == FUNCTION_VARIANCE) {   
    for (int q = 0; q < nsqpts; q++){
      wq = pc->quadrature(q, zq, yq);
      ffmean += wq*pc->basis(0,zq)*fvals[0*nsqpts+q]*fvals[0*nsqpts+q];
    }
  }
  delete [] zq;
  delete [] yq;
  if (moment_type == FUNCTION_MEAN) {
    return fmean;
  } else {
    return ffmean - fmean*fmean;
  }
}

void TACSKSStochasticFunction::getElementSVSens( int elemIndex, TACSElement *element,
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

        // Get the element basis class
        TACSElementBasis *basis = delem->getElementBasis();

        if (basis){

          for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){

            double pt[3];
            double weight = basis->getQuadraturePoint(i, pt);
      
            TacsScalar quantity = 0.0;
            delem->evalPointQuantity(elemIndex,
                                     this->quantityType,
                                     time, i, pt,
                                     Xpts, uq, udq, uddq,
                                     &quantity);

            TacsScalar Xd[9], J[9];
            TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);
          
            TacsScalar ksPtWeight = exp(ksWeight*(quantity - maxValue[j*nsqpts+q]))/ksSum[j*nsqpts+q];
            ksPtWeight *= weight*detJ*wt;
            if (moment_type == FUNCTION_VARIANCE){
              ksPtWeight *= 2.0*fvals[j*nsqpts+q];
            }
            TacsScalar dfdq = ksPtWeight;
            delem->addPointQuantitySVSens(elemIndex,
                                          this->quantityType,
                                          time,
                                          alpha, beta, gamma,
                                          i, pt,
                                          Xpts, uq, udq, uddq,
                                          &dfdq, dfduj);

          } // spatial integration

        }

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

void TACSKSStochasticFunction::addElementDVSens( int elemIndex, TACSElement *element,
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
  
  for (int j = 0; j < 1; j++){

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
        TACSElementBasis *basis = delem->getElementBasis();
        
        if (basis){

          for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){

            double pt[3];
            double weight = basis->getQuadraturePoint(i, pt);

            TacsScalar quantity = 0.0;
            delem->evalPointQuantity(elemIndex,
                                     this->quantityType,
                                     time, i, pt,
                                     Xpts, uq, udq, uddq,
                                     &quantity);        
          
            TacsScalar Xd[9], J[9];
            TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

            TacsScalar dfdq = wt*weight*detJ*exp(ksWeight*(quantity - maxValue[j*nsqpts+q]))/ksSum[j*nsqpts+q];
            if (moment_type == 1){
              dfdq *= 2.0*fvals[j*nsqpts+q];
            }
            delem->addPointQuantityDVSens( elemIndex, 
                                           this->quantityType,
                                           time, scale,
                                           i, pt,
                                           Xpts, uq, udq, uddq, &dfdq, 
                                           dvLen, dfdxj ); 

          } // spatial integration
        
        }

      }
      
    } // end yloop

    // for (int n = 0; n < nnodes; n++){
    //   int lptr = n*dvpernode;
    //   int gptr = n*dvpernode + j*dvpernode;
    //   for (int i = 0; i < dvpernode; i++){
    //     dfdx[gptr+n] += dfdxj[lptr+n];
    //   }            
    // }

    // printf("check nodewise placement of derivatives");
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

