!> Monte Carlo swarm data generator. Writes transport data (like mobility, mean energy etc)
!! to the file output/[simulation name]_TD.txt
!!
!! Code written by Jannis Teunissen
!! Last update: 16-04-2014

program particle_swarm
   use m_config
   use m_particle_swarm
   use m_io

   implicit none

   integer, parameter :: dp = kind(0.0d0)
   character(len=80) :: swarm_name, cfg_name, gas_name, out_file

   logical           :: accurate
   integer           :: mm
   integer           :: info_cntr
   integer           :: i_fld, n_flds
   integer           :: n_swarms_min, n_swarms, swarm_size

   real(dp)          :: tau_coll, tau_swarm, growth_rate
   real(dp)          :: min_fld, max_fld, weight

   real(dp)          :: init_eV, scale_factor
   real(dp)          :: fld, factor, temp

   integer, parameter :: n_pc = 4
   type(PC_t) :: pcs(n_pc)
   real(dp), allocatable :: all_td(:, :)

   call initialize_all()

   ! Initialize variables
   info_cntr     = 0
   n_swarms_min = CFG_get_int("swarm_min_number")
   swarm_size = CFG_get_int("swarm_size")
   scale_factor = CFG_get_real("swarm_fld_scale_factor")
   init_eV = CFG_get_real("init_eV")
   n_flds = CFG_get_int("swarm_n_flds")
   min_fld = CFG_get_real("swarm_fld_range", 1)
   max_fld = CFG_get_real("swarm_fld_range", 2)
   call CFG_get("td_abs_acc", abs_acc)
   call CFG_get("td_rel_acc", rel_acc)

   allocate(all_td(n_flds, PM_num_td))
   call IO_write_coeffs("output/" // trim(swarm_name) // "_coeffs.txt")

   do i_fld = 1, n_flds
      factor   = dble(i_fld-1) / max(1, n_flds-1)
      temp     = 1.0_dp / scale_factor
      temp     = factor * (max_fld**temp - min_fld**temp) + min_fld**temp
      fld      = temp ** scale_factor

      call PM_get_swarm_data(pcs(1), fld, all_td(i_fld, :))

      write(*,*) ""
      do mm = 1, PM_num_td
         write(*, "(A25,E10.3,A,E8.2)") "   " // PM_td_names(mm), all_td(i_fld, mm), " +- ", &
              td_dev(mm)
      end do
   end do

   gas_name = CFG_get_string("gas_mixture_name")
   out_file = "output/" // trim(swarm_name) // "_td.txt"
   call IO_write_td(all_td, PM_td_names, trim(gas_name), out_file, 1, 2)
   print *, "Output written to " // trim(out_file)
   out_file = "output/" // trim(swarm_name) // "_td_cols.txt"
   call IO_write_td_cols(all_td, PM_td_names, trim(gas_name), out_file)
   print *, "Output written to " // trim(out_file)

contains

   !> Initializes everything needed for the simulation
   subroutine initialize_all()
      use m_gas
      use m_cross_sec

      integer                              :: nn, n_gas_comp
      character(len=100)              :: cs_file
      character(len=20), allocatable :: gas_comp_names(:)
      real(dp), allocatable                :: gas_comp_fracs(:)
      type(CS_type), allocatable           :: cross_secs(:)

      call create_sim_config()      ! Create default parameters for the simulation (routine contained below)
      call CFG_sort()               ! Sort parameters in config module

      do nn = 1, command_argument_count()
         call get_command_argument(nn, cfg_name)
         call CFG_read_file(trim(cfg_name))
      end do

      call CFG_get("swarm_name", swarm_name)
      call CFG_write("output/" // trim(swarm_name) // "_config.txt")

      n_gas_comp = CFG_get_size("gas_component_names")
      if (n_gas_comp /= CFG_get_size("gas_component_fractions")) &
           print *, "gas_component_names and gas_component_fracs have unequal size"
      allocate(gas_comp_names(n_gas_comp))
      allocate(gas_comp_fracs(n_gas_comp))

      call CFG_get("gas_component_names", gas_comp_names)
      call CFG_get("gas_component_fractions", gas_comp_fracs)

      ! Initialize gas and electric field module
      call GAS_initialize(gas_comp_names, gas_comp_fracs, &
           CFG_get_real("gas_pressure"), CFG_get_real("gas_temperature"))

      do nn = 1, n_gas_comp
         cs_file = CFG_get_string("gas_crosssec_files", nn)
         call CS_add_from_file("input/" // trim(cs_file), &
              trim(gas_comp_names(nn)), 1.0_dp, &
              gas_comp_fracs(nn) * GAS_get_number_dens(), &
              CFG_get_real("part_max_energy_eV"), cross_secs)
      end do

      call CS_write_summary("output/" // trim(swarm_name) // "_cs_summary.txt")

      call pcs(1)%initialize(UC_elec_mass, cross_secs, &
              CFG_get_int("part_lkptbl_size"), &
              CFG_get_real("part_max_energy_eV"), &
              CFG_get_int("part_max_number_of"))

      do nn = 2, n_pc
         pcs(nn) = pcs(1)
      end do

   end subroutine initialize_all

   !> Create the parameters and their default values for the simulation
   subroutine create_sim_config()

      ! General simulation parameters
      call CFG_add("swarm_fld_range", (/1.0e6_dp, 1.0e7_dp/), "The minimum required electric field")
      call CFG_add("swarm_n_flds", 10, "The minimum required electric field")
      call CFG_add("swarm_fld_scale_factor", 1.0_dp, "The scaling of the electric fields")
      call CFG_add("swarm_name", "my_sim", "The name of the swarm simulation")
      call CFG_add("swarm_min_number", 10, "The minimum number of swarms")
      call CFG_add("swarm_size", 1000, "The initial size of a swarm")

      call CFG_add("td_abs_acc", (/1.0_dp, 0.0_dp, 0.0_dp, 1.0e-2_dp, 1.0e1_dp, 1.0e1_dp, 0.0_dp, 0.0_dp/), &
           "The required absolute accuracies of the transport data")
      call CFG_add("td_rel_acc", (/1.0_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-1_dp, 1e-1_dp/), &
           "The required relative accuracies of the transport data")
      call CFG_add("init_eV", 1.0_dp, "The initial energy of particles")

      ! Gas parameters
      call CFG_add("gas_pressure", 1.0_DP, "The gas pressure (bar)")
      call CFG_add("gas_temperature", 293.0_DP, "The gas temperature (Kelvin)")
      call CFG_add("gas_mixture_name", "N2", "The name of the gas mixture used")
      call CFG_add("gas_component_names", (/"N2"/), "The names of the gases used in the simulation", .true.)
      call CFG_add("gas_crosssec_files", (/"cross_sections_nitrogen.txt"/), &
           & "The files in which to find cross section data for each gas", .true.)
      call CFG_add("gas_component_fractions", (/1.0_dp /), &
           & "The partial pressure of the gases (as if they were ideal gases)", .true.)
      ! call CFG_add("gas_crosssec_scaling", 1.0_DP, "Scale factor for the cross sections in the input data", .true.)

      ! Particle model related parameters
      call CFG_add("part_lkptbl_size", 20*1000, "The size of the lookup table for the collision rates")
      call CFG_add("part_max_number_of", 25*1000*1000, "The maximum number of particles allowed per task")
      call CFG_add("part_max_energy_eV", 1000.0_DP, "The maximum energy in eV for particles in the simulation")
   end subroutine create_sim_config

end program particle_swarm
