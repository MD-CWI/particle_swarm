program particle_swarm
  use m_config
  use m_particle_core
  use m_particle_swarm
  use m_io

  implicit none

  integer, parameter  :: dp = kind(0.0d0)
  type(CFG_t)         :: cfg
  character(len=80)   :: swarm_name
  integer             :: swarm_size
  type(PC_t)          :: pc    ! Particle data
  type(SWARM_field_t) :: field ! The field configuration
  type(SWARM_td_t)    :: td(SWARM_num_td)
  logical             :: dry_run, visualize_only

  call initialize_all(cfg)

  if (.not. dry_run) then
     call get_field_configuration(cfg, field)
     call CFG_get(cfg, "swarm_size", swarm_size)
     call SWARM_initialize(pc, cfg, td, field)

     if (visualize_only) then
        call SWARM_visualize(pc, swarm_size, cfg)
     else
        call SWARM_get_data(pc, swarm_size, td)
        call SWARM_print_results(td)
     end if
  end if

contains

  !> Initializes everything needed for the simulation
  subroutine initialize_all(cfg)
    use m_gas
    use m_cross_sec
    use m_units_constants
    use omp_lib
    type(CFG_t), intent(inout)     :: cfg
    integer                        :: nn, tbl_size, max_num_part
    integer                        :: swarm_size, n_gas_comp, n_gas_frac
    integer                        :: rng_seed(4), num_threads
    real(dp)                       :: pressure, temperature, max_ev
    real(dp)                       :: magnetic_field, electric_field, tmp
    character(len=200)             :: cs_file, output_dir, tmp_name
    character(len=20)              :: particle_mover
    character(len=20), allocatable :: gas_names(:)
    real(dp), allocatable          :: gas_fracs(:)
    type(CS_t), allocatable        :: cross_secs(:)
    logical                        :: consecutive_run

    ! Create default parameters for the simulation (routine contained below)
    call create_sim_config(cfg)
    call CFG_sort(cfg)

    call CFG_update_from_arguments(cfg)
    call CFG_get(cfg, "swarm_name", swarm_name)

    call CFG_get(cfg, "num_threads", num_threads)
    if (num_threads > 0) then
       call omp_set_num_threads(num_threads)
    end if

    call CFG_get(cfg, "dry_run", dry_run)
    call CFG_get(cfg, "visualize_only", visualize_only)

    call CFG_get_size(cfg, "gas_components", n_gas_comp)
    call CFG_get_size(cfg, "gas_fractions", n_gas_frac)
    if (n_gas_comp /= n_gas_frac) &
         error stop "gas_components and gas_component_fracs have unequal size"
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
       call CFG_get(cfg, "particle_max_energy_ev", max_ev)

       do nn = 1, n_gas_comp
          call CS_add_from_file(trim(cs_file), &
               trim(gas_names(nn)), gas_fracs(nn) * &
               GAS_number_dens, max_ev, cross_secs)
       end do

       call CS_write_summary(cross_secs, &
            trim(output_dir) // "/" // trim(swarm_name) // "_cs_summary.txt")

       call CFG_get(cfg, "particle_lkptbl_size", tbl_size)
       call CFG_get(cfg, "swarm_size", swarm_size)

       call CFG_get(cfg, "particle_rng_seed", rng_seed)
       if (all(rng_seed == 0)) rng_seed = get_random_seed()

       if (visualize_only) then
          call CFG_get(cfg, "visualize_max_particles", max_num_part)
       else
          ! Allocate storage for 8 times the swarm size. There are checks in place
          ! to make sure it cannot grow to such a large size.
          max_num_part = 8 * swarm_size
       end if

       call pc%initialize(UC_elec_mass, max_num_part, rng_seed)
       call pc%use_cross_secs(max_ev, tbl_size, cross_secs)

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
            trim(tmp_name) // "_lt.dat", rng_seed)
    end if

    call CFG_get(cfg, "particle_mover", particle_mover)
    call CFG_get(cfg, "magnetic_field", magnetic_field)
    call CFG_get(cfg, "electric_field", electric_field)

    select case (trim(particle_mover))
    case ("analytic")
       if (electric_field > 10.0_dp * magnetic_field * UC_lightspeed) then
          ! Negligible B-field, use simplified approximation
          pc%particle_mover => SWARM_particle_mover_simple
       else
          pc%particle_mover => SWARM_particle_mover_analytic
       end if
    case ("boris")
       pc%particle_mover => SWARM_particle_mover_boris
       ! Limit time step to boris_dt_factor / cyclotron freq.
       call CFG_get(cfg, "boris_dt_factor", tmp)
       pc%dt_max = tmp * 2 * UC_pi / abs(magnetic_field * UC_elec_q_over_m)
    case ("verlet")
       pc%particle_mover => PC_verlet_advance
    case default
       error stop "Incorrect particle mover selected"
    end select

  end subroutine initialize_all

  !> Create the parameters and their default values for the simulation
  subroutine create_sim_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "num_threads", -1, &
         "Number of OpenMP threads to use (default: all)")
    call CFG_add(cfg, "consecutive_run", .false., &
         "True means: use data from a previous run with the same name")
    call CFG_add(cfg, "dry_run", .false., &
         "True means: only write simulation settings, no results")
    call CFG_add(cfg, "visualize_only", .false., &
         "True means: only generate files to visualize swarm")
    call CFG_add(cfg, "visualize_max_particles", 100*1000, &
         "Maximum number of particles for swarm visualization")
    call CFG_add(cfg, "visualize_rotate_Ez", .false., &
         "Rotate results so that E points in the z-direction")
    call CFG_add(cfg, "visualize_init_v0", [0.0_dp, 0.0_dp, 0.0_dp], &
         "Initial velocity of particles (if relax_swarm is false)")
    call CFG_add(cfg, "visualize_relax_swarm", .true., &
         "If true, relax the swarm first")

    ! Settings for visualization
    call CFG_add(cfg, "visualize_end_time", 1e-9_dp, &
         "End time for visualizing swarm")
    call CFG_add(cfg, "visualize_dt_output", 5e-11_dp, &
         "Time between writing visualization files")
    call CFG_add(cfg, "visualize_base_name", "particle_data", &
         "Base file name for visualization files")

    ! General simulation parameters
    call CFG_add(cfg, "electric_field", 1.0e7_dp, &
         "The electric field")
    call CFG_add(cfg, "magnetic_field", 0.0_dp, &
         "The electric field")
    call CFG_add(cfg, "field_angle_degrees", 90.0_dp, &
         "The angle between the electric and magnetic field")
    call CFG_add(cfg, "swarm_name", "my_sim", &
         "The name of the swarm simulation")
    call CFG_add(cfg, "swarm_size", 1000, &
         "The initial size of a swarm")
    call CFG_add(cfg, "output_dir", "output", &
         "The output directory (include no trailing slash!)")

    call CFG_add(cfg, "acc_velocity_sq", [5.0e-3_dp, 0.0_dp], &
         "The required rel/abs accuracy of the velocity squared")
    call CFG_add(cfg, "acc_velocity", [5.0e-3_dp, 0.0_dp], &
         "The required rel/abs accuracy of the velocity")
    call CFG_add(cfg, "acc_diffusion", [1.0e-2_dp, 0.0_dp], &
         "The required rel/abs accuracy of the diffusion coeff.")
    call CFG_add(cfg, "acc_alpha", [5.0e-3_dp, 0.0_dp], &
         "The required rel/abs accuracy of the ionization coeff.")

    ! Gas parameters
    call CFG_add(cfg, "gas_pressure", 1.0_dp, &
         "The gas pressure (bar)")
    call CFG_add(cfg, "gas_temperature", 300.0_dp, &
         "The gas temperature (Kelvin)")
    call CFG_add(cfg, "gas_components", (/"N2"/), &
         "The names of the gases used in the simulation", .true.)
    call CFG_add(cfg, "gas_file", "input/cross_sections_nitrogen.txt", &
         "The file in which to find cross section data")
    call CFG_add(cfg, "gas_fractions", (/1.0_dp /), &
         & "The partial pressure of the gases (as if they were ideal gases)", .true.)

    ! Particle model related parameters
    call CFG_add(cfg, "particle_rng_seed", [0, 0, 0, 0], &
         "Seed for random numbers. If all zero generate seed from clock.")
    call CFG_add(cfg, "particle_lkptbl_size", 10000, &
         "The size of the lookup table for the collision rates")
    call CFG_add(cfg, "particle_max_energy_ev", 500.0_dp, &
         "The maximum energy in eV for particles in the simulation")
    call CFG_add(cfg, "particle_mover", "verlet", &
         "Which particle mover to use. Options: analytic, verlet, boris")
    call CFG_add(cfg, "boris_dt_factor", 0.1_dp, &
         "The maximum time step in terms of the cylotron frequency")
  end subroutine create_sim_config

  function get_random_seed() result(seed)
    integer :: seed(4)
    integer :: time, i

    call system_clock(time)
    do i = 1, 4
       seed(i) = ishftc(time, i*8)
    end do
  end function get_random_seed

  subroutine get_field_configuration(cfg, field)
    use m_units_constants
    type(CFG_t), intent(inout)         :: cfg
    type(SWARM_field_t), intent(inout) :: field
    real(dp)                           :: electric_field, degrees

    ! Set electric and magnetic fields
    call CFG_get(cfg, "magnetic_field", field%Bz)
    call CFG_get(cfg, "electric_field", electric_field)
    call CFG_get(cfg, "field_angle_degrees", degrees)

    field%angle_deg = degrees
    field%Ez        = electric_field * cos(UC_pi * degrees / 180)
    field%Ey        = electric_field * sin(UC_pi * degrees / 180)

    field%B_vec = [0.0_dp, 0.0_dp, field%Bz]
    field%E_vec = [0.0_dp, field%Ey, field%Ez]

    field%omega_c   = -UC_elec_charge * field%Bz / UC_elec_mass
    field%ExB_drift = [field%Ey * field%Bz, 0.0_dp, 0.0_dp] / &
         max(epsilon(1.0_dp), field%Bz**2)

  end subroutine get_field_configuration

end program particle_swarm
