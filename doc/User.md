# User documentation

This documentation explains what is inside of the file User.py. Hence, all variables of the program is here. The variables are sorted in different dictionnaries.

## Units used

The lenght is in µm.

The time is in second.

The force is in µN.

The  mass is in kg.

## dict_algorithm

#### dt_PF
It is the time step used in phase-field simulations.
#### n_t_PF
It is approximatively the number of iteration during the phase-field simulation (an implicit scheme is used).
#### factor_distribution_etai
description in coming...
#### MovePF_selector
It is the method used to move the phase-fields. The following methods are available : DeconstructRebuild.
#### dt_DEM
It is the time step used in discrete element model simulations. It can be computed from a critical time step.
#### factor_neighborhood
description in coming...
#### Spring_type
It is the method used to compute the interaction grain-grain and grain-wall. The following methods are availablle : Ponctual.
#### i_DEM_stop
It is a stop criteria. During one discrete element model simulation, the number of iterations can not reach over i_DEM_stop.
#### Ecin_ratio
It is a stop criteria. Description in coming...
#### n_window_stop
It is a stop criteria. It works with dk0_stop and dy_box_max_stop. This variable defines the number of iterations considered in the extract window.
#### dk0_stop
It is a stop criteria. It works with n_window_stop. A window of n_window_stop iterations is extracted. The maximum and the miminum of the k0 in this window is compared. If the amplitude is lower than dk0_stop, the convergence criteria is partially verified. To obtain the total convergence criteria, the condition obtained with dy_box_max_stop must be verified.
#### dy_box_max_stop
It is a stop criteria. It works with n_window_stop. A window of n_window_stop iterations is extracted. The maximum and the miminum of the upper wall position in this window is compared. If the amplitude is lower than dy_box_max_stop, the convergence criteria is partially verified. To obtain the total convergence criteria, the condition obtained with dk0_stop must be verified.
#### n_t_PFDEM
It is the number of iterations done during the PFDEM simulation. The DEM simulation and the PF simulation are called n_t_PFDEM times.
#### np_proc
It is the number of processor used to run this simulation.
#### Debug
It is a boolean to plot configuration before and after the DEM simulation.
#### Debug_DEM
It is a boolean to plot configuration during the DEM simulation.
#### i_print_plot
It is the frequency of the print and the plot (plot only if Debug_DEM is True) during the DEM simulation.
#### SaveData
It is a boolean to save data.
#### main_folder_name
It is the name of the main folder where data is saved (if SaveData is True).
#### template_simulation_name
It is the name of the folder inside the main folder where data is saved (if SaveData is True). An index i (=1,2,...) is added at the end of the template.

## dict_geometry
#### N_grain_disk
It is the number of grains asked created on a disk basis.
#### L_R 
It is the different radius of grain distribution. This list must be sorted from the largest value to the smallest. This list must be sorted in agreement with L_percentage_R. The grain distribution is applied to the grain created on a disk basis.
#### L_percentage_R
It is the distribution of the different radius. This list must be sorted in agreement with L_R. The total sum must be 1. The grain distribution is applied to the grain created on a disk basis.
#### N_grain_square
It is the number of grains asked created on a square basis.
#### L_Dimension
It is the different dimension of grain distribution. This list must be sorted from the largest value to the smallest. This list must be sorted in agreement with L_percentage_Dimension. The grain distribution is applied to the grain created on a square basis. The dimension is the lenght of the square.
#### L_percentage_Dimension
It is the distribution of the different radius. This list must be sorted in agreement with L_Dimension. The total sum must be 1. The grain distribution is applied to the grain created on a square basis.

## dict_ic

#### n_generation
It is the number of generation used to create a initial condition. A generation is based on grain creation then loading until settlement.
#### factor_ymax_box
description in coming...
#### N_test_max
It is the maximum of number of tries to create a grain without an overlap with already created grains. The center of the grain are selected randomly.
#### i_DEM_stop_IC
It is a stop criteria. During discrete element model simulation tom compute initial configuration, the number of iterations can not reach over i_DEM_stop_IC.
#### Debug_DEM_IC
It is a boolean to plot configuration during the DEM simulation of the initial configuration.
#### i_print_plot_IC
It is the frequency of the print and the plot (plot only if Debug_DEM_IC is True) during the DEM simulation of the initial configuration.
#### dt_DEM_IC
It is the time step used in discrete element model simulation for the initial configuration. It can be computed from a critical time step.
#### Ecin_ratio_IC
It is a stop criteria. Description in coming...
#### factor_neighborhood_IC
description in coming...
#### i_update_neighborhoods_gen
description in coming...
#### i_update_neighborhoods_com
description in coming...

## dict_material

#### Y
It is the Young Modulus of the grain material.
#### nu
It is the Poisson's ratio of the grain material.
#### rho
It is the density of the grain material.
#### mu_friction_gg
It is the friction coefficient used in a Coulomb law for the interaction grain-grain.
#### mu_friction_gw
It is the friction coefficient used in a Coulomb law for the interaction wall-grain.
#### M 
It is the mobility parameter used in phase-field equations.
#### kc_pf
It is the gradient coefficient used in phase-field equations.
#### w
It is the width of the interface gradient used in phase-field simulations.
#### double_well_height
It is the energy barrier used in the double well formulation of the phase-field simulations. 

## dict_sample

#### x_box_min
It is the position of the wall on the left of the sample.
#### x_box_max
It is the position of the wall on the right of the sample.
#### y_box_min
It is the position of the lower wall of the sample.
#### x_min
It is the left domain limit of the phase field simulation.
#### x_max
It is the right domain limit of the phase field simulation.
#### n_x
It is the number of nodes in the x-direction.
#### y_min
It is the lower domain limit of the phase field simulation.
#### y_max
It is the upper domain limit of the phase field simulation.
#### n_y
It is the number of nodes in the y-direction.


## dict_sollicitations

#### Vertical_confinement_Force
It is the force applied on the upper wall to confine the sample.
#### gravity
It is the gravity applied on the sample. Can be 0 µm/s².
#### frac_dissolved
It is the fraction of particle dissolved. Must be between (or equal) 0 and 1.
#### Dissolution_Energy
It is the energy added in phase-field equations to dissolve the grains.
## dict_tracker
