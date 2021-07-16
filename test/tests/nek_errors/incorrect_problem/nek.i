[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 4
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [dummy]
    type = Diffusion
    variable = u
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  exodus = true
[]

[Postprocessors]
  [temp_average]
    type = NekVolumeAverage
    field = temperature
  []
[]
