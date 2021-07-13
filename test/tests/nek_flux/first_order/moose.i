[Mesh]
  type = FileMesh
  file = ../../postprocessors/meshes/pyramid.exo

  # MOOSE auxkernels to compute heat flux are limited to MONOMIAL bases
  # due to how material properties are stored. For comparison with Nek,
  # just use a much finer mesh
  uniform_refine = 3
[]

[Variables]
  [dummy]
  []
[]

[Kernels]
  [d]
    type = Diffusion
    variable = dummy
  []
[]

[AuxVariables]
  [temperature]
  []
  [flux]
    order = CONSTANT
    family = MONOMIAL
  []
  [difference]
  []
[]

[ICs]
  [temperature]
    type = FunctionIC
    variable = temperature
    function = temperature
  []
[]

[AuxKernels]
  [flux]
    type = DiffusionFluxAux
    variable = flux
    diffusion_variable = temperature
    component = normal
    diffusivity = thermal_conductivity
    boundary = '1 2 3 4 5 6 7 8'
    check_boundary_restricted = false
  []
[]

[Materials]
  [k]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '3.5'
  []
[]

[Executioner]
  type = Steady
[]

[Functions]
  [temperature]
    type = ParsedFunction
    value = 'exp(x)+sin(y)+x*y*z'
  []
[]

[Postprocessors]
  [flux_integral]
    type = SideIntegralVariablePostprocessor
    variable = flux
    boundary = '1 2 3 4 5 6 7 8'
  []
  [max_flux]
    type = ElementExtremeValue
    variable = flux
    value_type = max
  []
  [min_flux]
    type = ElementExtremeValue
    variable = flux
    value_type = min
  []
[]

[Outputs]
  exodus = true
  hide = 'dummy'
[]
