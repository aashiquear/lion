#
# This input file an initial test of the phase field wetting model
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 30
  ny = 15
  nz = 0
  xmax = 10
  ymax = 5
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 4
[]

[Variables]
  # We are solving the Cahn-Hilliard equation with the split form
  # The split form requires two variables
  [./phi]
    # The phase variable, 1 = liquid, -1 = vapor
  [../]
  [./chem_pot]
    # The chemical potential, mu
  [../]
[]

[AuxVariables]
  # These variables define the solid and the wetting behavior
  [./a]
    # Solid variable, 1 = liquid/vaport, -1 = solid
  [../]
  [./A]
    # Wetting parameter, > 0 is hydrophilic, < 0 is hydrophobic
    initial_condition = 0.9
  [../]
[]

[Functions]
  [./rectangle_image]
    type = ImageFunction
    upper_value = -1
    lower_value = 1
    file = rectangles.png
    threshold = 300
  [../]
[]

[ICs]
  active = 'alphaIC phi_IC'
  [./alphaIC]
    # This initial condition will be constant through the simulation
    # 1=fluid,gas; -1=solid
    x1 = 0
    x2 = 10
    y1 = 0
    y2 = 0.025
    inside = -1
    outside = 1
    variable = a
    type = BoundingBoxIC
  [../]
  [./alpha_rough_IC]
    type = FunctionIC
    function = rectangle_image
    variable = a
  [../]
  [./phi_IC]
    variable = phi
    type = SmoothCircleIC
    invalue = 1
    outvalue = -1
    radius = 2.2
    x1 = 5
    y1 = 2.7
    int_width = 0.1
  [../]
  [./phi_IC2]
    variable = phi
    type = SmoothCircleIC
    invalue = 1
    outvalue = -1
    radius = 3.2
    x1 = 5
    y1 = -2.4
    int_width = 0.1
  [../]
[]

[Kernels]
  # Kernels define different pieces of our equation we are solving
  # dot(phi) = laplacian(-gamma*laplacian(phi) + dV/dphi)
  # but we solve it with the split form, so we have two equations
  # variable = mu: dot(phi) = laplacian(mu)
  # variable = phi: mu = -gamma*laplacian(phi) + dV/dphi
  [./c_dot]
    # dot(phi)
    type = CoupledTimeDerivative
    variable = chem_pot
    v = phi
  [../]
  [./c_res]
    # mu = -gamma*laplacian(phi) + dV/dphi
    type = SplitCHParsed
    variable = phi
    f_name = V
    kappa_name = gamma
    w = chem_pot
  [../]
  [./w_res]
    # laplacian(mu)
    type = SplitCHWRes
    variable = chem_pot
    mob_name = M # M isn't actually used in this model so we set it to 1
    args = phi
  [../]
[]

[Materials]
  active = 'V gamma M_const'
  [./V]
    # Calculates dV/dphi
    type = DerivativeParsedMaterial
    block = 0
    constant_expressions = 2
    function = '0.5*(1 + a)*0.25*(1 - phi^2)^2 + 0.5*(1 - a)*(K/2)*(phi - A)^2'
    outputs = exodus
    args = 'phi A a'
    constant_names = K
    f_name = V
    derivative_order = 2
  [../]
  [./gamma]
    type = GenericConstantMaterial
    block = 0
    prop_names = gamma
    prop_values = 0.000325
  [../]
  [./M_const]
    type = GenericConstantMaterial
    prop_names = 'M'
    prop_values = '1'
  [../]
  [./M]
    type = DerivativeParsedMaterial
    block = 0
    args = 'a phi'
    constant_names = 'M0 lb'
    constant_expressions = '1 0.0'
    outputs = exodus
    f_name = M
    derivative_order = 2
    function = 'ub:=if(a>0.9,1,0.05); if(phi<-1, lb*M0, if(phi>1, ub*M0, (M0*((lb + ub) + phi*(ub - lb))/2)))'
  [../]
[]

[BCs]
[]

[Preconditioning]
  [./off_diag_coupling]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./total_phi]
    type = ElementIntegralVariablePostprocessor
    variable = phi
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm         31   preonly   lu      1'
  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 20
  nl_rel_tol = 1e-8
  end_time = 1000.0
  nl_abs_tol = 1e-9
  dtmax = 5
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.75
    dt = 0.001
    growth_factor = 1.2
    iteration_window = 2
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    refine_fraction = 0.9
    coarsen_fraction = 0.05
    max_h_level = 4
    weight_names = 'phi chem_pot'
    weight_values = '1 0.1'
  [../]
[]

[Outputs]
  execute_on = 'initial timestep_end'
  exodus = true
  print_perf_log = true
  interval = 5
  file_base = test
[]
