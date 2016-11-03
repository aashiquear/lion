/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifndef GRADIENTDIRECTION_H
#define GRADIENTDIRECTION_H

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

// Forward Declaration
class GradientDirection;


template<>
InputParameters validParams<GradientDirection>();

/**
 * Compute convection term for phase field simulations
 */
class GradientDirection : public DerivativeMaterialInterface<JvarMapInterfaceBase<Kernel> >
{
public:
  GradientDirection(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _property;

  const MaterialProperty<Real> & _dpropertydu;

  /// number of coupled variables
  const unsigned int _nvar;

  /// @{ Mobility derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dpropertydarg;

};

#endif // GRADIENTDIRECTION_H
