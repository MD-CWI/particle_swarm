program particle_swarm
  use omp_lib
  use m_config
  use m_particle_core
  use m_particle_swarm
  use m_io

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  type(CFG_t)        :: cfg
  character(len=80)  :: swarm_name
  integer            :: mm
  integer            :: n_swarms_min, swarm_size
  real(dp)           :: fld
  type(PC_t)         :: pc               ! Particle data
  type(SWARM_acc_t)  :: acc              ! Accuracy reqs.
  real(dp)           :: td(SWARM_num_td) ! Transport data
  real(dp)           :: td_dev(SWARM_num_td) ! Transport data error
  logical            :: dry_run

  call initialize_all(cfg)
  call CFG_get(cfg, "dry_run", dry_run)

  if (.not. dry_run) then
     call CFG_get(cfg, "swarm_min_number", n_swarms_min)
     call CFG_get(cfg, "swarm_size", swarm_size)
     call CFG_get(cfg, "swarm_fld", fld)
     call CFG_get(cfg, "acc_energy", acc%energy)
     call CFG_get(cfg, "acc_mobility", acc%mobility)
     call CFG_get(cfg, "acc_diffusion", acc%diffusion)

     call SWARM_get_data(pc, fld, swarm_size, &
          n_swarms_min, acc, td, td_dev)

     do mm = 1, SWARM_num_td
        write(*, "(A25,E10.3,E9.2)") "   " // &
             SWARM_td_names(mm), td(mm), td_dev(mm)
     end do
  end if

contains

  !> Initializes everything needed for the simulation
  subroutine initialize_all(cfg)
    use m_gas
    use m_cross_sec
    use m_units_constants
    type(CFG_t), intent(inout)      :: cfg
    integer                         :: nn, tbl_size, max_num_part
    integer                         :: swarm_size, n_gas_comp, n_gas_frac
    real(dp)                        :: pressure, temperature, max_ev
    character(len=200)              :: cs_file, output_dir
    character(LEN=200)              :: cfg_name, prev_name, tmp_name
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

    ! Initialize gas and electric field module
    call GAS_initialize(gas_names, gas_fracs, pressure, temperature)

    call CFG_get(cfg, "consecutive_run", consecutive_run)
    call CFG_get(cfg, "output_dir", output_dir)

    if (.not. consecutive_run) then
       tmp_name = trim(output_dir) // "/" // trim(swarm_name) // "_config.txt"
       print *, "Writing configuration to ", trim(tmp_name)
       call CFG_write(cfg, trim(tmp_name))

       call CFG_get(cfg, "gas_file", cs_file)
       call CFG_get(cfg, "part_max_energy_ev", max_ev)

       do nn = 1, n_gas_comp
          call CS_add_from_file(trim(cs_file), &
               trim(gas_names(nn)), gas_fracs(nn) * &
               GAS_number_dens, max_ev, cross_secs)
       end do

       call CS_write_summary(cross_secs, &
            trim(output_dir) // "/" // trim(swarm_name) // "_cs_summary.txt")

       print *, "Initializing particle model", 1
       call CFG_get(cfg, "part_lkptbl_size", tbl_size)
       call CFG_get(cfg, "swarm_size", swarm_size)

       ! Allocate storage for 4 times the swarm size. There are actually checks
       ! in place to make sure it cannot grow to such a large size.
       max_num_part = 4 * swarm_size

       call pc%initialize(UC_elec_mass, cross_secs, &
            tbl_size, max_ev, max_num_part, get_random_seed())

       print *, "--------------------"
       print *, "Gas information"
       write(*, fmt="(A10,A3,E9.3)") "Temp. (K)", " - ", temperature
       print *, ""
       print *, "Component - fraction"
       do nn = 1, n_gas_comp
          write(*, fmt="(A10,A3,E9.3)") trim(gas_names(nn)), " - ", &
               gas_fracs(nn)
       end do
       print *, ""

       tmp_name = trim(output_dir) // "/" // trim(swarm_name) // "_coeffs.txt"
       print *, "Writing colrate table (as text) to ", trim(tmp_name)
       call IO_write_coeffs(pc, trim(tmp_name))

       tmp_name = trim(output_dir) // "/" // trim(swarm_name)
       print *, "Writing particle params (raw) to ", &
            trim(tmp_name) // "_params.dat"
       print *, "Writing particle lookup table (raw) to ", &
            trim(tmp_name) // "_lt.dat"
       call pc%to_file(trim(tmp_name) // "_params.dat", &
            trim(tmp_name) // "_lt.dat")
       print *, "--------------------"

    else ! Restarted run (can only change field!)
       tmp_name = trim(output_dir) // "/" // trim(swarm_name)
       call pc%init_from_file(trim(tmp_name) // "_params.dat", &
            trim(tmp_name) // "_lt.dat", get_random_seed())
    end if

  end subroutine initialize_all

  !> Create the parameters and their default values for the simulation
  subroutine create_sim_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "consecutive_run", .false., &
         "True means: use data from a previous run with the same name")
    call CFG_add(cfg, "dry_run", .false., &
         "True means: only write simulation settings, no results'")

    ! General simulation parameters
    call CFG_add(cfg, "swarm_fld", 1.0e7_dp, &
         "The electric field")
    call CFG_add(cfg, "swarm_name", "my_sim", &
         "The name of the swarm simulation")
    call CFG_add(cfg, "swarm_min_number", 10, &
         "The minimum number of swarms")
    call CFG_add(cfg, "swarm_size", 1000, &
         "The initial size of a swarm")
    call CFG_add(cfg, "output_dir", "output", &
         "The output directory (include no trailing slash!)")

    call CFG_add(cfg, "acc_energy", [1.0e-2_dp, 0.0_dp], &
         "The required rel/abs accuracy of the energy")
    call CFG_add(cfg, "acc_mobility", [1.0e-2_dp, 0.0_dp], &
         "The required rel/abs accuracy of the mobility")
    call CFG_add(cfg, "acc_diffusion", [1.0e-2_dp, 0.0_dp], &
         "The required rel/abs accuracy of the diffusion coeff.")

    ! Gas parameters
    call CFG_add(cfg, "gas_pressure", 1.0_dp, &
         "The gas pressure (bar)")
    call CFG_add(cfg, "gas_temperature", 300.0_dp, &
         "The gas temperature (Kelvin)")
    call CFG_add(cfg, "gas_mixture_name", "N2", &
         "The name of the gas mixture used")
    call CFG_add(cfg, "gas_components", (/"N2"/), &
         "The names of the gases used in the simulation", .true.)
    call CFG_add(cfg, "gas_file", "input/cross_sections_nitrogen.txt", &
         "The file in which to find cross section data")
    call CFG_add(cfg, "gas_fractions", (/1.0_dp /), &
         & "The partial pressure of the gases (as if they were ideal gases)", .true.)

    ! Particle model related parameters
    call CFG_add(cfg, "part_lkptbl_size", 10000, &
         "The size of the lookup table for the collision rates")
    call CFG_add(cfg, "part_max_energy_ev", 500.0_dp, &
         "The maximum energy in eV for particles in the simulation")
  end subroutine create_sim_config

  function get_random_seed() result(seed)
    integer :: seed(4)
    integer :: time, i

    call system_clock(time)
    do i = 1, 4
       seed(i) = ishftc(time, i*8)
    end do
  end function get_random_seed

end program particle_swarm
