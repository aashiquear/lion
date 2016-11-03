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
#include "MatTimeDerivative.h"

template<>
InputParameters validParams<MatTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addRequiredParam<MaterialPropertyName>("mat_prop", "Name of the property (scalar) to multiply the MatTimeDerivative kernel with");
  params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  return params;
}

MatTimeDerivative::MatTimeDerivative(const InputParameters & parameters) :
    DerivativeMaterialInterface<JvarMapInterfaceBase<TimeDerivative > >(parameters),
    _coeff(getMaterialProperty<Real>("mat_prop")),
    _dcoeffdu(getMaterialPropertyDerivative<Real>("mat_prop", _var.name())),
    _nvar(_coupled_moose_vars.size()),
    _dcoeffdarg(_nvar)
{
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    MooseVariable *ivar = _coupled_moose_vars[i];
    _dcoeffdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", ivar->name());
  }
}

Real
MatTimeDerivative::computeQpResidual()
{
  return _coeff[_qp] * TimeDerivative::computeQpResidual();
}

Real
MatTimeDerivative::computeQpJacobian()
{
  return _coeff[_qp] * TimeDerivative::computeQpJacobian() +
  _dcoeffdu[_qp] * _phi[_j][_qp] * TimeDerivative::computeQpResidual();
}

Real
MatTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  // get the coupled variable jvar is referring to
  unsigned int cvar;

  if (!mapJvarToCvar(jvar, cvar))
    return 0.0;

  return (*_dcoeffdarg[cvar])[_qp] * _phi[_j][_qp] * TimeDerivative::computeQpResidual();
}
