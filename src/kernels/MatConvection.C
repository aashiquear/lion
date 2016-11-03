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
#include "MatConvection.h"

template<>
InputParameters validParams<MatConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<MaterialPropertyName>("mat_prop", "Name of the property (scalar) to multiply the MatConvection kernel with");
  params.addRequiredParam<RealVectorValue>("driving_vector", "Driving vector");
  params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  return params;
}

MatConvection::MatConvection(const InputParameters & parameters) :
    DerivativeMaterialInterface<JvarMapInterfaceBase<Kernel> >(parameters),
    _conv_prop(getMaterialProperty<Real>("mat_prop")),
    _dconv_propdu(getMaterialPropertyDerivative<Real>("mat_prop", _var.name())),
    _driving_vector(getParam<RealVectorValue>("driving_vector")),
    _nvar(_coupled_moose_vars.size()),
    _dconv_propdarg(_nvar)
{
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    MooseVariable *ivar = _coupled_moose_vars[i];
    _dconv_propdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", ivar->name());
  }
}

Real
MatConvection::computeQpResidual()
{
  return _test[_i][_qp] * _conv_prop[_qp] * _driving_vector * _grad_u[_qp];
}

Real
MatConvection::computeQpJacobian()
{
  return _test[_i][_qp] * _driving_vector * _conv_prop[_qp] * _grad_phi[_j][_qp] +
  _test[_i][_qp] * _dconv_propdu[_qp] * _driving_vector * _grad_u[_qp] * _phi[_j][_qp];
}

Real
MatConvection::computeQpOffDiagJacobian(unsigned int jvar)
{
  // get the coupled variable jvar is referring to
  unsigned int cvar;

  if (!mapJvarToCvar(jvar, cvar))
    return 0.0;

  return _test[_i][_qp] * _driving_vector * _phi[_j][_qp] * (*_dconv_propdarg[cvar])[_qp] * _grad_u[_qp];
}
