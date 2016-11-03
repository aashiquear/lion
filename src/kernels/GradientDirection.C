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
#include "GradientDirection.h"

template<>
InputParameters validParams<GradientDirection>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<MaterialPropertyName>("mat_prop", "Name of the property (scalar) to multiply the GradientDirection kernel with");
  params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  return params;
}

GradientDirection::GradientDirection(const InputParameters & parameters) :
    DerivativeMaterialInterface<JvarMapInterfaceBase<Kernel> >(parameters),
    _property(getMaterialProperty<Real>("mat_prop")),
    _dpropertydu(getMaterialPropertyDerivative<Real>("mat_prop", _var.name())),
    _nvar(_coupled_moose_vars.size()),
    _dpropertydarg(_nvar)
{
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    MooseVariable *ivar = _coupled_moose_vars[i];
    _dpropertydarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", ivar->name());
  }
}

Real
GradientDirection::computeQpResidual()
{
  Real res = 0;
  if (_grad_u[_qp].norm() > 0)
    res = _property[_qp] * _grad_u[_qp] / _grad_u[_qp].norm() * _grad_test[_i][_qp];

  return res;
}

Real
GradientDirection::computeQpJacobian()
{
  Real jac = 0;
  if (_grad_u[_qp].norm() > 0)
  {
    RealGradient dgradunorm = _grad_u[_qp]/_grad_u[_qp].norm();

    jac = _property[_qp] * _grad_phi[_j][_qp] / _grad_u[_qp].norm() * _grad_test[_i][_qp]
    - _property[_qp] * _grad_u[_qp] / _grad_u[_qp].norm_sq() * dgradunorm * _grad_phi[_j][_qp] * _grad_test[_i][_qp]
    + _dpropertydu[_qp] * _phi[_j][_qp] * _grad_u[_qp] / _grad_u[_qp].norm() * _grad_test[_i][_qp];
  }

  return jac;
}

Real
GradientDirection::computeQpOffDiagJacobian(unsigned int jvar)
{
  // get the coupled variable jvar is referring to
  unsigned int cvar;

  if (!mapJvarToCvar(jvar, cvar))
    return 0.0;

  return (*_dpropertydarg[cvar])[_qp] * _phi[_j][_qp] * _grad_u[_qp] / _grad_u[_qp].norm() * _grad_test[_i][_qp];
}
