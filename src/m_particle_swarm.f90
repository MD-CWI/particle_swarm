module m_particle_swarm
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  !> The axis of gyration for electrons points in the z-direction (just as the
  !> magnetic field)
  real(dp), parameter :: SWARM_omega_unitvec(3) = [0.0_dp, 0.0_dp, 1.0_dp]

  ! These are the basic swarm parameters that we measure, but some others that
  ! can be derived from these (e.g., the mean energy or mobility) are printed as
  ! well.
  integer, parameter :: SWARM_num_td  = 10
  integer, parameter :: ix_alpha      = 1
  integer, parameter :: ix_eta        = 2
  integer, parameter :: ix_ionization = 3
  integer, parameter :: ix_attachment = 4
  integer, parameter :: ix_coll_rate  = 5
  integer, parameter :: ix_flux_v     = 6
  integer, parameter :: ix_bulk_v     = 7
  integer, parameter :: ix_flux_diff  = 8
  integer, parameter :: ix_bulk_diff  = 9
  integer, parameter :: ix_flux_v2    = 10


  !> Indices of ionization collisions
  integer, allocatable :: ionization_colls(:)

  !> Indices of attachment collisions
  integer, allocatable :: attachment_colls(:)

  !> Type for storing transport data
  type SWARM_td_t
     character(len=20)     :: description    !< Description
     character(len=20)     :: unit           !< Unit
     integer               :: n_dim          !< Dimension of td(:)
     integer               :: n_measurements !< Number of measurements
     real(dp), allocatable :: val(:)         !< Actual values
     real(dp), allocatable :: var(:)         !< Estimated variance (sum)
     real(dp)              :: rel_acc        !< Req. relative accuracy
     real(dp)              :: abs_acc        !< Req. absolute accuracy
  end type SWARM_td_t

  !> Type storing the field configuration for a swarm
  type SWARM_field_t
     real(dp) :: B_vec(3)     !< Magnetic field vector
     real(dp) :: E_vec(3)     !< Electric field vector
     real(dp) :: angle_deg    !< Angle between E and B in degrees
     real(dp) :: Bz           !< Magnetic field along z-axis (T)
     real(dp) :: Ez           !< Electric field along z-axis (V/m)
     real(dp) :: Ey           !< Electric field along y-axis (V/m)
     real(dp) :: omega_c      !< Gyration frequency
     real(dp) :: ExB_drift(3) !< ExB drift velocity
  end type SWARM_field_t

  !> Type to collect particle statistics / properties
  type part_stats_t
     type(SWARM_field_t) :: field       !< The field configuration
     integer             :: n_samples   !< Number of samples
     real(dp)            :: coll_rate   !< Collision rate
     real(dp)            :: i_rate      !< Ionization rate
     real(dp)            :: a_rate      !< Attachment rate
     real(dp)            :: flux_v(3)   !< Mean velocity
     real(dp)            :: bulk_v(3)   !< Mean bulk velocity
     real(dp)            :: bulk_dif(3) !< Bulk diffusion
     real(dp)            :: flux_v2(3)  !< Mean square velocity
     real(dp)            :: cov_xv(3)   !< Covariance x,v
     real(dp)            :: x_prev(3)   !< Previous value of <x>
     real(dp)            :: x2_prev(3)  !< Previous value of <x^2>
  end type part_stats_t

  type(SWARM_field_t), protected :: SWARM_field

  public :: SWARM_td_t
  public :: SWARM_field_t
  public :: SWARM_num_td

  ! Public routines from this module
  public :: SWARM_initialize
  public :: SWARM_visualize
  public :: SWARM_get_data
  public :: SWARM_print_results
  public :: SWARM_particle_mover_analytic
  public :: SWARM_particle_mover_simple
  public :: SWARM_particle_mover_boris

contains

  subroutine swarm_initialize(pc, cfg, tds, field)
    use m_config
    use m_cross_sec
    type(PC_t), intent(in)          :: pc
    type(CFG_t), intent(inout)      :: cfg
    type(SWARM_td_t), intent(inout) :: tds(:)
    type(SWARM_field_t), intent(in) :: field
    real(dp)                        :: rel_abs_acc(2)
    integer                         :: n, i_i, i_a

    SWARM_field = field

    call init_td(tds(ix_flux_v2), 3, "flux_v2", "(m/s)^2")
    call init_td(tds(ix_flux_v), 3, "flux_v", "m/s")
    call init_td(tds(ix_bulk_v), 3, "bulk_v", "m/s")
    call init_td(tds(ix_flux_diff), 3, "flux_D", "m^2/s")
    call init_td(tds(ix_bulk_diff), 3, "bulk_D", "m^2/s")
    call init_td(tds(ix_alpha), 1, "alpha", "1/m")
    call init_td(tds(ix_eta), 1, "eta", "1/m")
    call init_td(tds(ix_coll_rate), 1, "coll_rate", "1/s")
    call init_td(tds(ix_ionization), 1, "ionization", "1/s")
    call init_td(tds(ix_attachment), 1, "attachment", "1/s")

    ! Get accuracy requirements
    call CFG_get(cfg, "acc_velocity_sq", rel_abs_acc)
    tds(ix_flux_v2)%rel_acc = rel_abs_acc(1)
    tds(ix_flux_v2)%abs_acc = rel_abs_acc(2)

    call CFG_get(cfg, "acc_velocity", rel_abs_acc)
    tds(ix_flux_v)%rel_acc = rel_abs_acc(1)
    tds(ix_flux_v)%abs_acc = rel_abs_acc(2)
    tds(ix_bulk_v)%rel_acc = rel_abs_acc(1)
    tds(ix_bulk_v)%abs_acc = rel_abs_acc(2)

    call CFG_get(cfg, "acc_diffusion", rel_abs_acc)
    tds(ix_flux_diff)%rel_acc = rel_abs_acc(1)
    tds(ix_flux_diff)%abs_acc = rel_abs_acc(2)
    tds(ix_bulk_diff)%rel_acc = rel_abs_acc(1)
    tds(ix_bulk_diff)%abs_acc = rel_abs_acc(2)

    call CFG_get(cfg, "acc_alpha", rel_abs_acc)
    tds(ix_alpha)%rel_acc = rel_abs_acc(1)
    tds(ix_alpha)%abs_acc = rel_abs_acc(2)

    ! Get indices of all ionization and attachment collisions
    n = count(pc%colls(1:pc%n_colls)%type == CS_ionize_t)
    allocate(ionization_colls(n))
    n = count(pc%colls(1:pc%n_colls)%type == CS_attach_t)
    allocate(attachment_colls(n))

    i_i = 0
    i_a = 0

    do n = 1, pc%n_colls
       if (pc%colls(n)%type == CS_ionize_t) then
          i_i = i_i + 1
          ionization_colls(i_i) = n
       else if (pc%colls(n)%type == CS_ionize_t) then
          i_a = i_a + 1
          attachment_colls(i_a) = n
       end if
    end do
  end subroutine swarm_initialize

  subroutine init_td(td, n_dim, description, unit)
    type(SWARM_td_t), intent(inout) :: td
    integer, intent(in)             :: n_dim
    character(len=*), intent(in)    :: description
    character(len=*), intent(in)    :: unit

    allocate(td%val(n_dim))
    allocate(td%var(n_dim))

    td%description    = description
    td%unit           = unit
    td%n_dim          = n_dim
    td%n_measurements = 0
    td%val            = 0.0_dp
    td%var            = 0.0_dp
    td%rel_acc        = 0.0_dp
    td%abs_acc        = huge(1.0_dp)
  end subroutine init_td

  !> Advance a swarm over time
  subroutine swarm_advance(pc, tau, desired_num_part, growth_rate)
    type(PC_t), intent(inout)        :: pc
    real(dp), intent(in)             :: tau
    integer, intent(in)              :: desired_num_part
    real(dp), intent(in)             :: growth_rate
    integer                          :: n, n_steps
    real(dp)                         :: dt

    ! Sometimes a swarm can rapidly grow or shink in time. Therefore we advance
    ! the swarm in steps, so that we can resize it if necessary.
    n_steps = ceiling(tau * growth_rate / log(2.0_dp))
    n_steps = max(n_steps, 1)
    dt      = tau/n_steps

    do n = 1, n_steps
       call pc%advance_openmp(dt)
       call resize_swarm(pc, desired_num_part)
    end do
  end subroutine swarm_advance

  !> Move swarm to that its center of mass is at the origin
  subroutine recenter_swarm(pc, avg_x)
    type(PC_t), intent(inout) :: pc
    real(dp), intent(out)     :: avg_x(3)
    real(dp)                  :: sum_x(3)

    call pc%compute_vector_sum(get_position, sum_x)
    avg_x = sum_x / pc%get_num_sim_part()

    call pc%translate(-avg_x)
  end subroutine recenter_swarm

  !> Place all particles at the origin
  subroutine collapse_swarm(pc)
    type(PC_t), intent(inout) :: pc
    call pc%loop_iopart(reset_part_x)
  end subroutine collapse_swarm

  subroutine get_position(part, x)
    type(PC_part_t), intent(in) :: part
    real(dp), intent(out)       :: x(:)
    x = part%x
  end subroutine get_position

  subroutine reset_part_x(part)
    type(PC_part_t), intent(inout) :: part
    part%x = 0.0_dp
  end subroutine reset_part_x

  subroutine resize_swarm(pc, new_size)
    use m_random
    type(PC_t), intent(inout) :: pc
    integer, intent(in) :: new_size

    integer :: ix, cur_size
    real(dp) :: chance
    type(PC_part_t) :: part

    cur_size = pc%get_num_sim_part()

    if (new_size >= 2 * cur_size) then
       do
          ! Double particles
          do ix = 1, cur_size
             part = pc%particles(ix)
             call pc%add_part(part)
          end do
          cur_size = cur_size * 2
          if (new_size < 2 * cur_size) exit
       end do
    else if (new_size < cur_size/2) then
       ! Reduce number of particles
       chance = new_size / real(cur_size, dp)
       do ix = 1, cur_size
          if (pc%rng%unif_01() > chance) call pc%remove_part(ix)
       end do
       call pc%clean_up()
    end if
  end subroutine resize_swarm

  subroutine initialize_particle_stats(ps, pc)
    type(part_stats_t), intent(inout) :: ps
    type(pc_t), intent(inout)         :: pc
    integer                           :: i, num_part

    ps%n_samples = 0
    ps%flux_v    = 0.0_dp
    ps%bulk_v    = 0.0_dp
    ps%bulk_dif  = 0.0_dp
    ps%flux_v2   = 0.0_dp
    ps%cov_xv    = 0.0_dp
    ps%coll_rate = 0.0_dp
    ps%i_rate    = 0.0_dp
    ps%a_rate    = 0.0_dp

    call recenter_swarm(pc, ps%x_prev)

    num_part = pc%get_num_sim_part()
    do i = 1, 3
       ps%x2_prev(i) = sum(pc%particles(1:num_part)%x(i)**2) / num_part
    end do
  end subroutine initialize_particle_stats

  subroutine update_particle_stats(pc, ps, dt)
    type(PC_t), intent(inout)         :: pc
    type(part_stats_t), intent(inout) :: ps
    real(dp), intent(in)              :: dt
    integer                           :: ix, num_part
    real(dp)                          :: inv_n_samples
    real(dp)                          :: v(3), v2(3), corr_fac, x2(3)
    real(dp)                          :: mean_x(3), bulk_v(3), bulk_dif(3)
    real(dp)                          :: coll_rates(PC_max_num_coll)
    type(PC_part_t)                   :: part

    num_part = pc%get_num_sim_part()
    corr_fac = num_part/(num_part-1.0_dp)

    call recenter_swarm(pc, mean_x)

    ! Update bulk velocity
    bulk_v    = mean_x/dt
    ps%bulk_v = ps%bulk_v + (bulk_v - ps%bulk_v) * &
         num_part / (ps%n_samples + num_part)

    ! Update bulk diffusion coefficient
    do ix = 1, 3
       x2(ix) = sum(pc%particles(1:num_part)%x(ix)**2) / num_part
    end do
    bulk_dif = (x2 - ps%x2_prev) / (2 * dt)

    ps%x2_prev = x2
    ps%bulk_dif = ps%bulk_dif + (bulk_dif - ps%bulk_dif) * &
         num_part / (ps%n_samples + num_part)

    do ix = 1, num_part
       part          = pc%particles(ix)
       ps%n_samples  = ps%n_samples + 1
       inv_n_samples = 1.0_dp / ps%n_samples
       v             = part%v
       v2            = v**2
       ps%flux_v     = ps%flux_v + (v - ps%flux_v) * inv_n_samples

       ! This placing is intentional:
       ! http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
       ps%cov_xv     = ps%cov_xv + (v - ps%flux_v) * part%x * corr_fac

       ! ps%cov_v2_xv  = ps%cov_v2_xv + v2 * (v - ps%v) * part%x
       ps%flux_v2         = ps%flux_v2 + (v2 - ps%flux_v2) * inv_n_samples

       call pc%get_coll_rates(sqrt(sum(v2)), coll_rates(1:pc%n_colls))
       ps%coll_rate = ps%coll_rate + (sum(coll_rates(1:pc%n_colls)) - &
            ps%coll_rate) * inv_n_samples
       ps%i_rate = ps%i_rate + (sum(coll_rates(ionization_colls)) - &
            ps%i_rate) * inv_n_samples
       ps%a_rate = ps%a_rate + (sum(coll_rates(attachment_colls)) - &
            ps%a_rate) * inv_n_samples
    end do

  end subroutine update_particle_stats

  subroutine update_td_from_ps(tds, ps)
    use m_units_constants
    type(SWARM_td_t), intent(inout) :: tds(:)
    type(part_stats_t), intent(in)  :: ps

    call update_td(tds(ix_flux_v2), ps%flux_v2)
    call update_td(tds(ix_flux_v), ps%flux_v)
    call update_td(tds(ix_bulk_v), ps%bulk_v)
    call update_td(tds(ix_bulk_diff), ps%bulk_dif)
    call update_td(tds(ix_flux_diff), ps%cov_xv / ps%n_samples)
    call update_td(tds(ix_alpha), [ps%i_rate / norm2(ps%flux_v)])
    call update_td(tds(ix_eta), [ps%a_rate / norm2(ps%flux_v)])
    call update_td(tds(ix_coll_rate), [ps%coll_rate])
    call update_td(tds(ix_ionization), [ps%i_rate])
    call update_td(tds(ix_attachment), [ps%a_rate])

    ! td(7)     = abs(ps%v2_v(3) / (fld * ps%v2))           ! Energy mobility
    ! td(8)     = ps%cov_v2_xv(1) / (ps%n_samples * ps%v2)  ! Long. energy diffusion
  end subroutine update_td_from_ps

  subroutine update_td(td, new_val)
    type(SWARM_td_t), intent(inout) :: td
    real(dp), intent(in)            :: new_val(td%n_dim)
    real(dp)                        :: delta(td%n_dim), fac

    td%n_measurements = td%n_measurements + 1
    fac               = 1.0_dp / td%n_measurements
    delta             = new_val - td%val

    ! Update mean
    td%val = td%val + fac * delta

    ! Update variance
    td%var = td%var + delta * (new_val - td%val)
  end subroutine update_td

  !> Initialize a new swarm
  subroutine new_swarm(pc, swarm_size, init_eV, init_v0)
    use m_units_constants
    use m_random
    type(PC_t), intent(inout)       :: pc
    integer, intent(in)             :: swarm_size
    real(dp), intent(in), optional  :: init_eV    !< Initial energy (in eV)
    real(dp), intent(in), optional  :: init_v0(3) !< Initial velocity (in m/s)
    integer                         :: ll
    real(dp)                        :: x(3), v(3), a(3), vel

    call pc%remove_particles()

    x   = 0
    a   = SWARM_field%E_vec * UC_elec_q_over_m

    if (present(init_eV)) then
       vel = sqrt(2 * init_eV * UC_elem_charge/UC_elec_mass)

       do ll = 1, swarm_size
          v   = pc%rng%sphere(vel)
          call pc%create_part(x, v, a, 1.0_dp, 0.0_dp)
       end do
    else if (present(init_v0)) then
       do ll = 1, swarm_size
          call pc%create_part(x, init_v0, a, 1.0_dp, 0.0_dp)
       end do
    end if
  end subroutine new_swarm

  !> Produce files to visualize a swarm of electrons
  subroutine SWARM_visualize(pc, swarm_size, cfg)
    use m_config
    type(PC_t), intent(inout)  :: pc
    integer, intent(in)        :: swarm_size
    type(CFG_t), intent(inout) :: cfg
    character(len=200)         :: base_name
    integer                    :: n, n_steps, verbose
    real(dp)                   :: t_end, dt_output, t
    real(dp)                   :: v0(3), rotmat(2,2), mean_x(3)
    logical                    :: rotate, relax_swarm
    type(part_stats_t)         :: ps

    call CFG_get(cfg, "visualize_end_time", t_end)
    call CFG_get(cfg, "visualize_dt_output", dt_output)
    call CFG_get(cfg, "visualize_base_name", base_name)
    call CFG_get(cfg, "visualize_rotate_Ez", rotate)
    call CFG_get(cfg, "visualize_init_v0", v0)
    call CFG_get(cfg, "visualize_relax_swarm", relax_swarm)
    call CFG_get(cfg, "verbose", verbose)

    n_steps = nint(t_end/dt_output)

    if (n_steps > 1000*1000) &
         error stop "SWARM_visualize: more than 1e6 steps"
    if (swarm_size > size(pc%particles)) &
         error stop "SWARM_visualize: increase SWARM_max_particles"

    ! Rotate around x-axis so E points in the z direction
    if (rotate) then
       t = -SWARM_field%angle_deg * acos(-1.0_dp) / 180_dp
       rotmat(:, 1) = [cos(t), -sin(t)]
       rotmat(:, 2) = [sin(t), cos(t)]
       v0(2:3)      = matmul(rotmat, v0(2:3))
    else
       rotmat(:, 1) = [1.0_dp, 0.0_dp]
       rotmat(:, 2) = [0.0_dp, 1.0_dp]
    end if

    if (relax_swarm) then
       call create_swarm(pc, ps, swarm_size, verbose)
       call recenter_swarm(pc, mean_x)
    else
       call new_swarm(pc, swarm_size, init_v0=v0)
    end if

    call write_particles(pc, trim(base_name), 0.0_dp, 0, rotmat)

    do n = 1, n_steps
       call pc%advance_openmp(dt_output)
       call write_particles(pc, trim(base_name), n * dt_output, n, rotmat)
    end do
  end subroutine SWARM_visualize

  subroutine write_particles(pc, base_name, time, cntr, rotmat)
    type(PC_t), intent(inout)    :: pc
    character(len=*), intent(in) :: base_name
    character(len=200)           :: file_name
    real(dp), intent(in)         :: time
    integer, intent(in)          :: cntr
    real(dp), intent(in)         :: rotmat(2,2)
    integer                      :: i
    integer, parameter           :: my_unit = 300
    real(dp), allocatable        :: pdata(:, :)
    real(dp), save               :: rmin(3), rmax(3)

    write(file_name, "(A,I6.6,A)") base_name // "_", cntr, ".txt"
    open(my_unit, file=trim(file_name))

    allocate(pdata(6, pc%n_part))
    pdata(1, :) = pc%particles(1:pc%n_part)%x(1)
    pdata(2, :) = pc%particles(1:pc%n_part)%x(2)
    pdata(3, :) = pc%particles(1:pc%n_part)%x(3)
    pdata(4, :) = pc%particles(1:pc%n_part)%v(1)
    pdata(5, :) = pc%particles(1:pc%n_part)%v(2)
    pdata(6, :) = pc%particles(1:pc%n_part)%v(3)

    pdata(2:3, :) = matmul(rotmat, pdata(2:3, :))
    pdata(5:6, :) = matmul(rotmat, pdata(5:6, :))

    if (cntr == 1) then
       rmin = minval(pdata(1:3, :), dim=2)
       rmax = maxval(pdata(1:3, :), dim=2)
    else
       rmin = min(rmin, minval(pdata(1:3, :), dim=2))
       rmax = max(rmax, maxval(pdata(1:3, :), dim=2))
    end if

    ! Write header
    write(my_unit, "(A)") "# File with particle coordinates: x(3) v(3)"
    write(my_unit, "(A,E12.4)") "# time =", time
    write(my_unit, "(A,3E12.4)") "# rmin (up to now) =", rmin
    write(my_unit, "(A,3E12.4)") "# rmax (up to now) =", rmax

    do i = 1, pc%n_part
       write(my_unit, *) pdata(:, i)
    end do

    close(my_unit)
    print *, "Particles written to: ", trim(file_name)

  end subroutine write_particles

  subroutine SWARM_get_data(pc, swarm_size, tds, verbose, max_cpu_time)
    use iso_fortran_env, only: error_unit
    use m_units_constants

    type(PC_t), intent(inout)       :: pc
    integer, intent(in)             :: swarm_size
    type(SWARM_td_t), intent(inout) :: tds(:)
    integer, intent(in)             :: verbose
    real(dp), intent(in)            :: max_cpu_time

    integer, parameter  :: n_swarms_min = 10
    integer, parameter  :: n_swarms_max = 10000
    real(dp), parameter :: fac          = 0.5 * UC_elec_mass/UC_elem_charge
    integer             :: n_measurements
    integer             :: n, n_swarms, imax(1)
    real(dp)            :: dt, t_relax, t_measure, growth_rate
    real(dp)            :: rel_error(SWARM_num_td), t0, t1
    type(part_stats_t)  :: ps

    call create_swarm(pc, ps, swarm_size, verbose)

    call cpu_time(t0)

    ! Loop over the swarms until converged
    do n_swarms = 1, n_swarms_max
       ! Estimate the time scale for energy relaxation, given by:
       ! energy / (d/dt energy)
       t_relax        = fac * sum(ps%flux_v2) / &
            abs(dot_product(SWARM_field%E_vec, ps%flux_v))
       t_measure      = t_relax
       ! Aim for one measurement per collision time
       n_measurements = nint(t_measure * ps%coll_rate)
       ! Limit range of n_measurements
       n_measurements = min(10000, max(10, n_measurements))
       dt             = t_measure / n_measurements
       growth_rate    = abs(ps%i_rate - ps%a_rate)

       ! Advance over energy relaxation time
       call collapse_swarm(pc)
       call swarm_advance(pc, t_relax, swarm_size, growth_rate)

       call initialize_particle_stats(ps, pc)

       do n = 1, n_measurements
          call swarm_advance(pc, dt, swarm_size, growth_rate)
          call update_particle_stats(pc, ps, dt)
       end do

       call update_td_from_ps(tds, ps)

       if (verbose > 1) call SWARM_print_results(tds, verbose)

       if (n_swarms >= n_swarms_min) then
          ! Check whether the results are accurate enough
          call get_accuracy(tds, rel_error)
          if (verbose > 0) then
             imax = maxloc(rel_error)
             write(*, '(I6,A35,F12.4)') n_swarms, 'max rel. error (' // &
                  trim(tds(imax(1))%description) // ')', rel_error(imax(1))
          end if

          call cpu_time(t1)
          if (all(rel_error < 1.0_dp) .or. t1 - t0 > max_cpu_time) exit
       end if
    end do

    if (n_swarms == n_swarms_max + 1) then
       write(error_unit, *) "No convergence in ", n_swarms_max, "iterations"
       error stop "No convergence"
    end if

  end subroutine SWARM_get_data

  !> Compute uncertainty in transport data relative to requirements
  subroutine get_accuracy(tds, rel_error)
    type(SWARM_td_t), intent(in) :: tds(:) !< The transport data
    real(dp), intent(out)        :: rel_error(SWARM_num_td)
    real(dp)                     :: stddev, mean, var
    integer                      :: i, n
    real(dp), parameter          :: eps = 1e-100_dp

    do i = 1, SWARM_num_td
       n      = tds(i)%n_measurements
       ! Below we use the norm for vector-based transport data, so that
       ! orientation of the field should not matter
       mean   = norm2(tds(i)%val)
       var    = maxval(tds(i)%var / (n-1)) / n ! Maximal variance
       stddev = sqrt(var)                      ! Standard deviation of the mean

       rel_error(i) = min(stddev / max(tds(i)%rel_acc * mean, eps), &
            stddev / max(tds(i)%abs_acc, eps))
    end do
  end subroutine get_accuracy

  ! Create a swarm that is relaxed to the electric field
  subroutine create_swarm(pc, ps, swarm_size, verbose)
    use m_units_constants
    !> Data type which stores particle model
    type(PC_t), intent(inout)       :: pc
    !> Recorded particle statistics
    type(part_stats_t), intent(out) :: ps
    !> Number of electrons in swarm
    integer, intent(in)             :: swarm_size
    integer, intent(in)             :: verbose

    integer, parameter    :: frame_size    = 100
    real(dp), parameter   :: en_eV         = 0.1_dp
    integer, parameter    :: max_its_relax = 500
    integer, parameter    :: min_its_relax = 5
    integer               :: i, cntr
    real(dp)              :: en_hist(frame_size), t_hist(frame_size)
    real(dp)              :: mean_en, correl, stddev, tau

    ! Create a new swarm with 1 eV electrons
    call new_swarm(pc, swarm_size, init_eV=1.0_dp)

    ! An electron accelerating from zero velocity gains en_eV in this time (in
    ! the absence of a magnetic field). This time step is only used to determine
    ! when the swarm is approximately relaxed to the background field.
    tau = sqrt(0.5_dp * en_eV * UC_elec_mass / UC_elem_charge) / &
         norm2(SWARM_field%E_vec)
    if (verbose > 1) print *, "dt for energy relaxation", tau

    ! Create linear table with unit variance and zero mean
    do i = 1, frame_size
       t_hist(i) = (i - 0.5_dp * (frame_size+1)) * &
            sqrt(12.0_dp / (frame_size + frame_size**2))
    end do

    call initialize_particle_stats(ps, pc)

    ! Determine when the mean energy is relaxed, so that we can use this swarm
    do cntr = 1, max_its_relax
       do i = 1, frame_size
          call pc%advance_openmp(tau)
          call update_particle_stats(pc, ps, tau)
          call resize_swarm(pc, swarm_size)
          en_hist(i) = pc%get_mean_energy() / UC_elec_volt
       end do

       mean_en = sum(en_hist) / frame_size
       stddev = sqrt(sum((en_hist - mean_en)**2) / (frame_size-1))

       ! Compute correlation between en_hist and a line
       correl = sum((en_hist - mean_en) * t_hist) / (frame_size * stddev)
       if (verbose > 1) print *, cntr, "energy correlation", correl
       ! If the correlation is sufficiently small, exit
       if (cntr > min_its_relax .and. abs(correl) < 0.25_dp) exit
    end do
  end subroutine create_swarm

  subroutine SWARM_print_results(tds, verbose)
    use m_units_constants
    type(SWARM_td_t), intent(in) :: tds(:) !< The transport data
    integer, intent(in)          :: verbose
    integer                      :: i, i_dim
    real(dp)                     :: fac, tmp, std
    real(dp)                     :: energy, mu, rel_error(SWARM_num_td)
    logical                      :: magnetic_field_used
    character(len=2)             :: dimnames(3) = ["_x", "_y", "_z"]

    i = tds(1)%n_measurements
    fac = sqrt(1.0_dp / (i * (i-1)))

    magnetic_field_used = (abs(SWARM_field%Bz) > 0.0_dp)

    if (verbose > 0) &
         write(*, "(A15,4A12)") "name", "value", "sigma", "convergence", "unit"

    write(*, "(A15,2E12.4,A24)") "E", norm2(SWARM_field%E_vec), 0.0_dp, "V/m"
    if (magnetic_field_used) then
       write(*, "(A15,2E12.4)") "B", SWARM_field%Bz, 0.0_dp
       write(*, "(A15,2E12.4)") "angle", SWARM_field%angle_deg, 0.0_dp
       write(*, "(A15,2E12.4)") "omega_c", SWARM_field%omega_c, 0.0_dp
    end if

    ! mean energy
    tmp    = 0.5_dp * UC_elec_mass / UC_elec_volt
    energy = tmp * sum(tds(ix_flux_v2)%val)
    std    = tmp * sqrt(sum(tds(ix_flux_v2)%var)) * fac
    write(*, "(A15,2E12.4,A24)") "energy", energy, std, "eV"

    if (magnetic_field_used) then
       ! mobility parallel to E
       mu = -dot_product(tds(ix_flux_v)%val,  SWARM_field%E_vec) &
            / max(epsilon(1.0_dp), sum(SWARM_field%E_vec**2))
       std = fac * sqrt(dot_product(tds(ix_flux_v)%var,  SWARM_field%E_vec)) &
            / max(epsilon(1.0_dp), norm2(SWARM_field%E_vec)**1.5_dp)
       write(*, "(A15,2E12.4,A24)") "mu_E", mu, std, "m^2/(Vs)"

       ! mobility parallel to B
       if (abs(SWARM_field%Ez) > sqrt(epsilon(1.0_dp))) then
          mu = -tds(ix_flux_v)%val(3) / SWARM_field%Ez
          std = fac * sqrt(tds(ix_flux_v)%var(3)) / abs(SWARM_field%Ez)
       else
          mu = 0
          std = 0
       end if
       write(*, "(A15,2E12.4,A24)") "mu_B", mu, std, "m^2/(Vs)"

       ! mobility perpendicular to B (y-velocity over Ey)
       if (abs(SWARM_field%Ey) > sqrt(epsilon(1.0_dp))) then
          mu = -tds(ix_flux_v)%val(2) / SWARM_field%Ey
          std = fac * sqrt(tds(ix_flux_v)%var(2)) / abs(SWARM_field%Ey)
       else
          mu = 0
          std = 0
       end if
       write(*, "(A15,2E12.4,A24)") "mu_xB", mu, std, "m^2/(Vs)"

       ! mobility in ExB-direction (x-velocity over Ey = E_perp)
       if (abs(SWARM_field%Ey) > sqrt(epsilon(1.0_dp)) .and. &
            abs(SWARM_field%Bz) > sqrt(epsilon(1.0_dp))) then
          mu = tds(ix_flux_v)%val(1) / SWARM_field%Ey
          std = fac * sqrt(tds(ix_flux_v)%var(1)) / abs(SWARM_field%Ey)
       else
          mu = 0
          std = 0
       end if
       write(*, "(A15,2E12.4,A24)") "mu_ExB", mu, std, "m^2/(Vs)"
    else
       ! mobility parallel to E
       mu = -dot_product(tds(ix_flux_v)%val,  SWARM_field%E_vec) &
            / max(epsilon(1.0_dp), sum(SWARM_field%E_vec**2))
       std = fac * sqrt(dot_product(tds(ix_flux_v)%var,  SWARM_field%E_vec)) &
            / max(epsilon(1.0_dp), norm2(SWARM_field%E_vec)**1.5_dp)
       write(*, "(A15,2E12.4,A24)") "mu_flux", mu, std, "m^2/(Vs)"

       mu = -dot_product(tds(ix_bulk_v)%val,  SWARM_field%E_vec) &
            / max(epsilon(1.0_dp), sum(SWARM_field%E_vec**2))
       std = fac * sqrt(dot_product(tds(ix_bulk_v)%var,  SWARM_field%E_vec)) &
            / max(epsilon(1.0_dp), norm2(SWARM_field%E_vec)**1.5_dp)
       write(*, "(A15,2E12.4,A24)") "mu_bulk", mu, std, "m^2/(Vs)"
    end if

    call get_accuracy(tds, rel_error)

    ! The other swarm parameters
    do i = 1, SWARM_num_td
       if (tds(i)%n_dim == 1) then
          write(*, "(A15,2E12.4,F12.4,A12)") &
               trim(tds(i)%description), &
               tds(i)%val(1), sqrt(tds(i)%var(1)) * fac, &
               rel_error(i), trim(tds(i)%unit)
       else

          do i_dim = 1, tds(i)%n_dim
             write(*, "(A15,2E12.4,F12.4,A12)") &
                  trim(tds(i)%description) // dimnames(i_dim),&
                  tds(i)%val(i_dim), &
                  sqrt(tds(i)%var(i_dim)) * fac, &
                  rel_error(i), trim(tds(i)%unit)
          end do
       end if
    end do

  end subroutine SWARM_print_results

  !> Advance the particle position and velocity over time dt taking into account
  !> a constant electric and magnetic field, using the analytic solution.
  subroutine SWARM_particle_mover_analytic(self, part, dt)
    use m_units_constants
    class(PC_t), intent(in)        :: self
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: dt
    real(dp)                       :: rc(3), rc_rot(3), theta

    ! Determine rotation angle
    theta = dt * SWARM_field%omega_c

    ! First the motion perpendicular to the magnetic field

    ! Subtract the plasma drift velocity
    part%v = part%v - SWARM_field%ExB_drift

    ! Determine rc, the vector pointing from the center of gyration
    rc = cross_product(part%v, SWARM_omega_unitvec) / SWARM_field%omega_c

    ! Rotate position and velocity
    part%v = rotate_around_axis(part%v, SWARM_omega_unitvec, theta)
    rc_rot = rotate_around_axis(rc, SWARM_omega_unitvec, theta)

    ! Update the position with the change in guiding center
    part%x = part%x + (rc_rot - rc) + SWARM_field%ExB_drift * dt

    ! Add back the plasma drift velocity
    part%v = part%v + SWARM_field%ExB_drift

    ! Now the motion parallel to the magnetic field
    part%v(3) = part%v(3) + 0.5_dp * dt * UC_elec_q_over_m * SWARM_field%Ez
    part%x(3) = part%x(3) + dt * part%v(3)
    part%v(3) = part%v(3) + 0.5_dp * dt * UC_elec_q_over_m * SWARM_field%Ez

    ! Update time left
    part%t_left = part%t_left - dt
  end subroutine SWARM_particle_mover_analytic

  !> Advance the particle position and velocity over time dt taking into account
  !> a constant electric and magnetic field, using the analytic solution.
  subroutine SWARM_particle_mover_simple(self, part, dt)
    use m_units_constants
    class(PC_t), intent(in)        :: self
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: dt
    real(dp)                       :: a(3)

    ! Use simple approximation based on Lorentz force
    a = UC_elec_q_over_m * (SWARM_field%E_vec + &
         cross_product(part%v, SWARM_field%B_vec))

    part%v = part%v + 0.5_dp * dt * a
    part%x = part%x + dt * part%v

    a = UC_elec_q_over_m * (SWARM_field%E_vec + &
         cross_product(part%v, SWARM_field%B_vec))
    part%v = part%v + 0.5_dp * dt * a

    ! Update time left
    part%t_left = part%t_left - dt
  end subroutine SWARM_particle_mover_simple

  !> Advance the particle position and velocity over time dt taking into account
  !> a constant electric and magnetic field, using Boris' method.
  subroutine SWARM_particle_mover_boris(self, part, dt)
    use m_units_constants
    class(PC_t), intent(in)        :: self
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: dt
    real(dp) :: t_vec(3), tmp(3)

    ! Push the particle over dt/2
    part%x = part%x + 0.5_dp * dt * part%v ! Use the previous velocity
    part%v = part%v + 0.5_dp * dt * UC_elec_q_over_m * SWARM_field%E_vec

    ! Rotate the velocity
    tmp    = 0.5_dp * dt * SWARM_field%B_vec * UC_elec_q_over_m
    t_vec  = 2 * tmp / (1.d0 + (norm2(tmp))**2)
    tmp    = part%v + cross_product(part%v, tmp)
    tmp    = cross_product(tmp, t_vec)
    part%v = part%v + tmp

    ! Push the particle over dt/2
    part%v = part%v + 0.5_dp * dt * UC_elec_q_over_m * SWARM_field%E_vec
    part%x = part%x + 0.5_dp * dt * part%v ! Use the new velocity

    ! Update time left
    part%t_left = part%t_left - dt
  end subroutine SWARM_particle_mover_boris

  !> Rotate a vector around an axis over an angle theta
  pure function rotate_around_axis(vec, axis, theta) result(rot)
    real(dp), intent(in) :: vec(3)  !< Vector to rotate
    real(dp), intent(in) :: axis(3) !< Axis of rotation
    real(dp), intent(in) :: theta   !< Angle of rotation
    real(dp)             :: rot(3)

    ! Use Rodrigues' rotation formula, see
    ! https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
    rot = vec * cos(theta) + &
         cross_product(axis, vec) * sin(theta) + &
         axis * dot_product(axis, vec) * (1 - cos(theta))
  end function rotate_around_axis

  !> Return the cross product of vectors a and b
  pure function cross_product(a, b) result(vec)
    real(dp), intent(in) :: a(3)
    real(dp), intent(in) :: b(3)
    real(dp)             :: vec(3)

    vec(1) = a(2) * b(3) - a(3) * b(2)
    vec(2) = a(3) * b(1) - a(1) * b(3)
    vec(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

end module m_particle_swarm
