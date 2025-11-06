#include "TACSMutableElement3D.h"

// Constructor
TACSMutableElement3D::TACSMutableElement3D( TACSElement *_elem ){
  this->element = _elem;
  this->element->incref();
}

// Destructor
TACSMutableElement3D::~TACSMutableElement3D(){
  this->element->decref();
}

void TACSMutableElement3D::setDensity( TacsScalar _rho ){
  (void) _rho;
}
