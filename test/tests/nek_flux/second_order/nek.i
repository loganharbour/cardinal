[Problem]
  type = NekRSProblem
  incoming_BC = temperature
[]

[Mesh]
  type = NekRSMesh
  order = SECOND
  boundary = '1 2 3 4 5 6'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Postprocessors]
  [total_flux] # should approximately match flux_integral
    type = ElementIntegralVariablePostprocessor
    variable = avg_flux
  []
  [max_flux]
    type = ElementExtremeValue
    variable = avg_flux
    value_type = max
  []
  [min_flux]
    type = ElementExtremeValue
    variable = avg_flux
    value_type = min
  []
[]

[Outputs]
  exodus = true
[]
