[Mesh]
  type = GeneratedMesh
  dim = 2
  nx =
  ny =
  nz = 0
  xmin =
  xmax =
  ymin =
  ymax =
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

#[GlobalParams]
#  outputs = exodus
#[]

[Variables]

[]

[Kernels]

[]

[Materials]
  [./pfmobility]
    type = GenericConstantMaterial
    prop_names  = 'L kappa_c'
    prop_values =
    block = 0
  [../]
  [./chemical_free_energy]
    type = DerivativeParsedMaterial
    block = 0
    f_name = Fc
    args =
    constant_names       = 'barr_height '
    constant_expressions =
    material_property_names = e_dissolution
    function = barr_height*(
    enable_jit = true
    derivative_order = 2
  [../]
  [./e_dissolution_map]
    type = GenericFunctionMaterial
    block = 0
    prop_names = e_dissolution
    prop_values = e_dissolution_txt
  [../]
[]

[Functions]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'implicit-euler'

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type  -sub_pc_type '
  petsc_options_value = 'asm       lu'

  l_max_its = 100
  nl_max_its = 20
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-3
  nl_abs_tol = 1.0e-3
  start_time = 0.0
  end_time =

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt =
  [../]
[]

[Outputs]
  exodus = true
  [other]
    type = VTK
  []
[]
