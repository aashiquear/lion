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
#ifndef MATTIMEDERIVATIVE_H
#define MATTIMEDERIVATIVE_H

#include "TimeDerivative.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

// Forward Declaration
class MatTimeDerivative;


template<>
InputParameters validParams<MatTimeDerivative>();

/**
 * Compute convection term for phase field simulations
 */
class MatTimeDerivative : public DerivativeMaterialInterface<JvarMapInterfaceBase<TimeDerivative> >
{
public:
  MatTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _coeff;

  const MaterialProperty<Real> & _dcoeffdu;

  /// number of coupled variables
  const unsigned int _nvar;

  /// @{ Mobility derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dcoeffdarg;

};

#endif // MATTIMEDERIVATIVE_H
