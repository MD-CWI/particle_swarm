program particle_swarm
  use omp_lib
  use m_config
  use m_particle_core
  use m_particle_swarm
  use m_io

  implicit none

  integer, parameter      :: dp = kind(0.0d0)
  type(CFG_t)             :: cfg
  character(len=80)       :: swarm_name, cfg_name, gas_name, out_file

  integer                 :: mm
  integer                 :: info_cntr
  integer                 :: i_fld, n_flds
  integer                 :: my_pc, n_pc, pc_cntr
  integer                 :: n_swarms_min, swarm_size
  real(dp)                :: min_fld, max_fld
  real(dp)                :: init_eV, scale_factor
  real(dp)                :: factor, temp
  real(dp)                :: fld_range(2)
  type(PC_t), allocatable :: pcs(:)
  real(dp), allocatable   :: all_td(:, :)
  real(dp), allocatable   :: fld_list(:)
  real(dp)                :: abs_acc(SWARM_num_td), rel_acc(SWARM_num_td)

  n_pc = omp_get_max_threads()

  print *, "Number of threads", n_pc
  allocate(pcs(n_pc))

  call initialize_all(cfg)

  ! Initialize variables
  info_cntr    = 0
  call CFG_get(cfg, "swarm_min_number", n_swarms_min)
  call CFG_get(cfg, "swarm_size", swarm_size)
  call CFG_get(cfg, "swarm_fld_scale_factor", scale_factor)
  call CFG_get(cfg, "init_eV", init_eV)
  call CFG_get(cfg, "swarm_n_flds", n_flds)

  call CFG_get(cfg, "swarm_fld_range", fld_range)
  min_fld = fld_range(1)
  max_fld = fld_range(2)
  call CFG_get(cfg, "td_abs_acc", abs_acc)
  call CFG_get(cfg, "td_rel_acc", rel_acc)

  allocate(all_td(SWARM_num_td, n_flds))
  allocate(fld_list(n_flds))

  print *, "Writing used coefficients to ",&
       "output/" // trim(swarm_name) // "_coeffs.txt"
  call IO_write_coeffs(pcs(1), "output/" // trim(swarm_name) // "_coeffs.txt")

  do i_fld = 1, n_flds
     factor          = dble(i_fld-1) / max(1, n_flds-1)
     temp            = 1.0_dp / scale_factor
     temp            = factor * (max_fld**temp - min_fld**temp) + min_fld**temp
     fld_list(i_fld) = temp ** scale_factor
  end do

  print *, "Starting swarm calculations"
  pc_cntr    = 0

  !$omp parallel private(my_pc)
  my_pc = omp_get_thread_num() + 1
  !$omp do schedule(dynamic, 1)
  do i_fld = 1, n_flds
     call SWARM_get_data(pcs(my_pc), fld_list(i_fld), swarm_size, &
          n_swarms_min, abs_acc, rel_acc, all_td(:, i_fld))

     write(*,*) ""
     do mm = 1, SWARM_num_td
        write(*, "(A25,E10.3,A,E8.2)") "   " // SWARM_td_names(mm), &
             all_td(mm, i_fld)
     end do
  end do
  !$omp enddo
  !$omp end parallel

  call CFG_get(cfg, "gas_mixture_name", gas_name)
  out_file = "output/" // trim(swarm_name) // "_td.txt"
  call IO_write_td(all_td, SWARM_td_names, trim(gas_name), out_file, 1, 2)
  print *, "Output written to " // trim(out_file)
  out_file = "output/" // trim(swarm_name) // "_td_cols.txt"
  call IO_write_td_cols(all_td, SWARM_td_names, trim(gas_name), out_file)
  print *, "Output written to " // trim(out_file)

contains

  !> Initializes everything needed for the simulation
  subroutine initialize_all(cfg)
    use m_gas
    use m_cross_sec
    use m_units_constants
    type(CFG_t), intent(inout)     :: cfg
    integer                        :: nn, tbl_size, max_num_part
    integer                        :: n_gas_comp, n_gas_frac
    real(dp)                       :: pressure, temperature, max_ev
    character(len=100), allocatable :: cs_files(:)
    character(len=20), allocatable :: gas_comp_names(:)
    real(dp), allocatable          :: gas_comp_fracs(:)
    type(CS_t), allocatable        :: cross_secs(:)

    call create_sim_config(cfg)      ! Create default parameters for the
    ! simulation (routine contained below)
    call CFG_sort(cfg)               ! Sort parameters in config module

    do nn = 1, command_argument_count()
       call get_command_argument(nn, cfg_name)
       call CFG_read_file(cfg, trim(cfg_name))
    end do

    call CFG_get(cfg, "swarm_name", swarm_name)
    call CFG_write(cfg, "output/" // trim(swarm_name) // "_config.txt")

    call CFG_get_size(cfg, "gas_components", n_gas_comp)
    call CFG_get_size(cfg, "gas_fractions", n_gas_frac)
    if (n_gas_comp /= n_gas_frac) &
         print *, "gas_components and gas_component_fracs have unequal size"
    allocate(gas_comp_names(n_gas_comp))
    allocate(gas_comp_fracs(n_gas_comp))

    call CFG_get(cfg, "gas_components", gas_comp_names)
    call CFG_get(cfg, "gas_fractions", gas_comp_fracs)
    call CFG_get(cfg, "gas_temperature", temperature)
    call CFG_get(cfg, "gas_pressure", pressure)

    print *, "--------------------"
    print *, "Gas information"
    write(*, fmt="(A10,A3,E9.3)") "Temp. (K)", " - ", temperature
    print *, ""
    print *, "Component - fraction"
    do nn = 1, n_gas_comp
       write(*, fmt="(A10,A3,E9.3)") trim(gas_comp_names(nn)), " - ", &
            gas_comp_fracs(nn)
    end do
    print *, "--------------------"

    ! Initialize gas and electric field module
    call GAS_initialize(gas_comp_names, gas_comp_fracs, pressure, temperature)
    allocate(cs_files(n_gas_comp))
    call CFG_get(cfg, "gas_files", cs_files)
    call CFG_get(cfg, "part_max_energy_ev", max_ev)

    do nn = 1, n_gas_comp
       call CS_add_from_file("input/" // trim(cs_files(nn)), &
            trim(gas_comp_names(nn)), gas_comp_fracs(nn) * &
            GAS_get_number_dens(), max_ev, cross_secs)
    end do

    call CS_write_summary(cross_secs, &
         "output/" // trim(swarm_name) // "_cs_summary.txt")

    print *, "Initializing particle model", 1
    call CFG_get(cfg, "part_lkptbl_size", tbl_size)
    call CFG_get(cfg, "part_max_number_of", max_num_part)

    call pcs(1)%initialize(UC_elec_mass, cross_secs, &
         tbl_size, max_ev, max_num_part)

    do nn = 2, n_pc
       print *, "Initializing particle model", nn
       pcs(nn) = pcs(1)
    end do

  end subroutine initialize_all

  !> Create the parameters and their default values for the simulation
  subroutine create_sim_config(cfg)
    type(CFG_t), intent(inout) :: cfg
    ! General simulation parameters
    call CFG_add(cfg, "swarm_fld_range", (/1.0e6_dp, 1.0e7_dp/), "The minimum required electric field")
    call CFG_add(cfg, "swarm_n_flds", 10, "The minimum required electric field")
    call CFG_add(cfg, "swarm_fld_scale_factor", 1.0_dp, "The scaling of the electric fields")
    call CFG_add(cfg, "swarm_name", "my_sim", "The name of the swarm simulation")
    call CFG_add(cfg, "swarm_min_number", 10, "The minimum number of swarms")
    call CFG_add(cfg, "swarm_size", 10*1000, "The initial size of a swarm")

    call CFG_add(cfg, "td_abs_acc", (/1.0_dp, 0.0_dp, 0.0_dp, 1.0e-2_dp, 1.0e1_dp, 1.0e1_dp, 0.0_dp, 0.0_dp/), &
         "The required absolute accuracies of the transport data")
    call CFG_add(cfg, "td_rel_acc", (/1.0_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-1_dp, 1e-1_dp/), &
         "The required relative accuracies of the transport data")
    call CFG_add(cfg, "init_eV", 1.0_dp, "The initial energy of particles")

    ! Gas parameters
    call CFG_add(cfg, "gas_pressure", 1.0_DP, "The gas pressure (bar)")
    call CFG_add(cfg, "gas_temperature", 293.0_DP, "The gas temperature (Kelvin)")
    call CFG_add(cfg, "gas_mixture_name", "N2", "The name of the gas mixture used")
    call CFG_add(cfg, "gas_components", (/"N2"/), "The names of the gases used in the simulation", .true.)
    call CFG_add(cfg, "gas_files", (/"cross_sections_nitrogen.txt"/), &
         & "The files in which to find cross section data for each gas", .true.)
    call CFG_add(cfg, "gas_fractions", (/1.0_dp /), &
         & "The partial pressure of the gases (as if they were ideal gases)", .true.)
    ! call CFG_add(cfg, "gas_crosssec_scaling", 1.0_DP, "Scale factor for the cross sections in the input data", .true.)

    ! Particle model related parameters
    call CFG_add(cfg, "part_lkptbl_size", 1*1000, "The size of the lookup table for the collision rates")
    call CFG_add(cfg, "part_max_number_of", 1000*1000, "The maximum number of particles allowed per task")
    call CFG_add(cfg, "part_max_energy_ev", 500.0_dp, "The maximum energy in eV for particles in the simulation")
  end subroutine create_sim_config

end program particle_swarm
