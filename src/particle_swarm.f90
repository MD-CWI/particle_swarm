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
  integer                 :: n_swarms_min, swarm_size
  real(dp)                :: fld
  real(dp)                :: init_eV, scale_factor
  type(PC_t) :: pc
  real(dp) :: td(SWARM_num_td)
  real(dp) :: abs_acc(SWARM_num_td)
  real(dp) :: rel_acc(SWARM_num_td)

  call initialize_all(cfg)

  ! Initialize variables
  call CFG_get(cfg, "swarm_min_number", n_swarms_min)
  call CFG_get(cfg, "swarm_size", swarm_size)
  call CFG_get(cfg, "init_eV", init_eV)
  call CFG_get(cfg, "swarm_fld", fld)
  call CFG_get(cfg, "td_abs_acc", abs_acc)
  call CFG_get(cfg, "td_rel_acc", rel_acc)

  print *, "Writing used coefficients to ",&
       "output/" // trim(swarm_name) // "_coeffs.txt"
  ! call IO_write_coeffs(pc, "output/" // trim(swarm_name) // "_coeffs.txt")

  print *, "Starting swarm calculations"
     call SWARM_get_data(pc, fld, swarm_size, &
          n_swarms_min, abs_acc, rel_acc, td)
     do mm = 1, SWARM_num_td
        write(*, "(A25,E10.3)") "   " // SWARM_td_names(mm), td(mm)
     end do

contains

  !> Initializes everything needed for the simulation
  subroutine initialize_all(cfg)
    use m_gas
    use m_cross_sec
    use m_units_constants
    type(CFG_t), intent(inout)      :: cfg
    integer                         :: nn, tbl_size, max_num_part
    integer                         :: n_gas_comp, n_gas_frac
    real(dp)                        :: pressure, temperature, max_ev
    character(len=100), allocatable :: cs_files(:)
    character(LEN=100)              :: filename, prev_name, tmp_name
    character(len=20), allocatable  :: gas_names(:)
    real(dp), allocatable           :: gas_fracs(:)
    type(CS_t), allocatable         :: cross_secs(:)
    logical                         :: consecutive_run

    ! Create default parameters for the simulation (routine contained below)
    call create_sim_config(cfg)
    call CFG_sort(cfg)

    swarm_name = "swarm"
    prev_name = ""
    do nn = 1, command_argument_count()
       call get_command_argument(nn, cfg_name)
       call CFG_read_file(cfg, trim(cfg_name))

       call CFG_get(cfg, "swarm_name", tmp_name)
       if (tmp_name /= "" .and. tmp_name /= prev_name) &
            swarm_name = trim(swarm_name) // "_" // trim(tmp_name)
       prev_name = tmp_name
    end do

    call CFG_get_size(cfg, "gas_components", n_gas_comp)
    call CFG_get_size(cfg, "gas_fractions", n_gas_frac)
    if (n_gas_comp /= n_gas_frac) &
         print *, "gas_components and gas_component_fracs have unequal size"
    allocate(gas_names(n_gas_comp))
    allocate(gas_fracs(n_gas_comp))

    call CFG_get(cfg, "gas_components", gas_names)
    call CFG_get(cfg, "gas_fractions", gas_fracs)
    call CFG_get(cfg, "gas_temperature", temperature)
    call CFG_get(cfg, "gas_pressure", pressure)

    print *, "--------------------"
    print *, "Gas information"
    write(*, fmt="(A10,A3,E9.3)") "Temp. (K)", " - ", temperature
    print *, ""
    print *, "Component - fraction"
    do nn = 1, n_gas_comp
       write(*, fmt="(A10,A3,E9.3)") trim(gas_names(nn)), " - ", &
            gas_fracs(nn)
    end do
    print *, "--------------------"

    ! Initialize gas and electric field module
    call GAS_initialize(gas_names, gas_fracs, pressure, temperature)

    call CFG_get(cfg, "consecutive_run", consecutive_run)

    if (.not. consecutive_run) then
       call CFG_write(cfg, "output/" // trim(swarm_name) // "_config.txt")

       allocate(cs_files(n_gas_comp))
       call CFG_get(cfg, "gas_files", cs_files)
       call CFG_get(cfg, "part_max_energy_ev", max_ev)

       do nn = 1, n_gas_comp
          call CS_add_from_file("input/" // trim(cs_files(nn)), &
               trim(gas_names(nn)), gas_fracs(nn) * &
               GAS_number_dens, max_ev, cross_secs)
       end do

       call CS_write_summary(cross_secs, &
            "output/" // trim(swarm_name) // "_cs_summary.txt")

       print *, "Initializing particle model", 1
       call CFG_get(cfg, "part_lkptbl_size", tbl_size)
       call CFG_get(cfg, "part_max_number_of", max_num_part)

       call pc%initialize(UC_elec_mass, cross_secs, &
            tbl_size, max_ev, max_num_part)
       call pc%to_file(trim(swarm_name) // "_params.dat", &
            trim(swarm_name) // "_lt.dat")

    else                        ! Restarted run
       call pc%init_from_file(trim(swarm_name) // "_params.dat", &
            trim(swarm_name) // "_lt.dat")
    end if

  end subroutine initialize_all

  !> Create the parameters and their default values for the simulation
  subroutine create_sim_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "consecutive_run", .false., &
         "True means: use data from a previous run with the same name")

    ! General simulation parameters
    call CFG_add(cfg, "swarm_fld", 1.0e7_dp, &
         "The electric field")
    call CFG_add(cfg, "swarm_name", "my_sim", &
         "The name of the swarm simulation")
    call CFG_add(cfg, "swarm_min_number", 10, &
         "The minimum number of swarms")
    call CFG_add(cfg, "swarm_size", 1000, &
         "The initial size of a swarm")

    call CFG_add(cfg, "td_abs_acc", &
         (/1.0_dp, 0.0_dp, 0.0_dp, 1.0e-2_dp, 1.0e1_dp, 1.0e1_dp, 0.0_dp, 0.0_dp/), &
         "The required absolute accuracies of the transport data")
    call CFG_add(cfg, "td_rel_acc", &
         (/1.0_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-2_dp, 1e-1_dp, 1e-1_dp/), &
         "The required relative accuracies of the transport data")
    call CFG_add(cfg, "init_eV", 1.0_dp, &
         "The initial energy of particles")

    ! Gas parameters
    call CFG_add(cfg, "gas_pressure", 1.0_dp, &
         "The gas pressure (bar)")
    call CFG_add(cfg, "gas_temperature", 293.0_dp, &
         "The gas temperature (Kelvin)")
    call CFG_add(cfg, "gas_mixture_name", "N2", &
         "The name of the gas mixture used")
    call CFG_add(cfg, "gas_components", (/"N2"/), &
         "The names of the gases used in the simulation", .true.)
    call CFG_add(cfg, "gas_files", &
         (/"cross_sections_nitrogen.txt"/), &
         & "The files in which to find cross section data for each gas", .true.)
    call CFG_add(cfg, "gas_fractions", (/1.0_dp /), &
         & "The partial pressure of the gases (as if they were ideal gases)", .true.)

    ! Particle model related parameters
    call CFG_add(cfg, "part_lkptbl_size", 20*1000, &
         "The size of the lookup table for the collision rates")
    call CFG_add(cfg, "part_max_number_of", 1000*1000, &
         "The maximum number of particles allowed per task")
    call CFG_add(cfg, "part_max_energy_ev", 900.0_dp, &
         "The maximum energy in eV for particles in the simulation")
  end subroutine create_sim_config

end program particle_swarm
