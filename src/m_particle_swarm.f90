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
  integer, parameter :: SWARM_num_td = 6
  integer, parameter :: ix_alpha     = 1
  integer, parameter :: ix_eta       = 2
  integer, parameter :: ix_coll_rate = 3
  integer, parameter :: ix_diffusion = 4
  integer, parameter :: ix_vel       = 5
  integer, parameter :: ix_vel_sq    = 6

  !> Type for storing transport data
  type SWARM_td_t
     character(len=20)     :: description    !< Description
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
     type(SWARM_field_t) :: field     !< The field configuration
     integer             :: n_samples !< Number of samples
     real(dp)            :: coll_rate !< Collision rate
     real(dp)            :: i_rate    !< Ionization rate
     real(dp)            :: a_rate    !< Attachment rate
     real(dp)            :: v(3)      !< Mean velocity
     real(dp)            :: v2(3)     !< Mean square velocity
     real(dp)            :: cov_xv(3) !< Covariance x,v
  end type part_stats_t

  type(SWARM_field_t), protected :: SWARM_field

  public :: SWARM_td_t
  public :: SWARM_field_t
  public :: SWARM_num_td

  ! Public routines from this module
  public :: SWARM_initialize
  public :: SWARM_get_data
  public :: SWARM_print_results
  public :: SWARM_particle_mover_analytic
  public :: SWARM_particle_mover_simple
  public :: SWARM_particle_mover_boris

contains

  subroutine swarm_initialize(cfg, tds, field)
    use m_config
    type(CFG_t), intent(in)         :: cfg
    type(SWARM_td_t), intent(inout) :: tds(:)
    type(SWARM_field_t), intent(in) :: field
    real(dp)                        :: rel_abs_acc(2)

    SWARM_field = field

    call init_td(tds(ix_vel_sq), 3, "vel_sq")
    call init_td(tds(ix_vel), 3, "vel")
    call init_td(tds(ix_diffusion), 3, "diff")
    call init_td(tds(ix_alpha), 1, "alpha")
    call init_td(tds(ix_eta), 1, "eta")
    call init_td(tds(ix_coll_rate), 1, "coll_rate")

    ! Get accuracy requirements
    call CFG_get(cfg, "acc_velocity_sq", rel_abs_acc)
    tds(ix_vel_sq)%rel_acc = rel_abs_acc(1)
    tds(ix_vel_sq)%abs_acc = rel_abs_acc(2)

    call CFG_get(cfg, "acc_velocity", rel_abs_acc)
    tds(ix_vel)%rel_acc = rel_abs_acc(1)
    tds(ix_vel)%abs_acc = rel_abs_acc(2)

    call CFG_get(cfg, "acc_diffusion", rel_abs_acc)
    tds(ix_diffusion)%rel_acc = rel_abs_acc(1)
    tds(ix_diffusion)%abs_acc = rel_abs_acc(2)

    call CFG_get(cfg, "acc_alpha", rel_abs_acc)
    tds(ix_alpha)%rel_acc = rel_abs_acc(1)
    tds(ix_alpha)%abs_acc = rel_abs_acc(2)

  end subroutine swarm_initialize

  subroutine init_td(td, n_dim, description)
    type(SWARM_td_t), intent(inout) :: td
    integer, intent(in)             :: n_dim
    character(len=*), intent(in)    :: description

    allocate(td%val(n_dim))
    allocate(td%var(n_dim))

    td%description    = description
    td%n_dim          = n_dim
    td%n_measurements = 0
    td%val            = 0.0_dp
    td%var            = 0.0_dp
    td%rel_acc        = 0.0_dp
    td%abs_acc        = huge(1.0_dp)
  end subroutine init_td

  !> Advance a swarm over time
  subroutine swarm_advance(pc, tau, desired_num_part, growth_rate)
    type(PC_t), intent(inout) :: pc
    real(dp), intent(in)      :: tau
    integer, intent(in)       :: desired_num_part
    real(dp), intent(in)      :: growth_rate
    integer                   :: n, n_steps
    real(dp)                  :: dt

    ! Sometimes a swarm can rapidly grow or shink in time. Therefore we advance
    ! the swarm in steps, so that we can resize it if necessary.
    n_steps = ceiling(tau * growth_rate / log(2.0_dp))
    n_steps = max(n_steps, 1)
    dt      = tau/n_steps

    do n = 1, n_steps
       call pc%advance(dt)
       call resize_swarm(pc, desired_num_part)
    end do
  end subroutine swarm_advance

  !> Move swarm to that its center of mass is at the origin
  subroutine recenter_swarm(pc)
    type(PC_t), intent(inout) :: pc
    real(dp)                  :: sum_x(3), avg_x(3)

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
             part = pc%get_part(ix)
             call pc%add_part(part)
          end do
          cur_size = cur_size * 2
          if (new_size < 2 * cur_size) exit
       end do
    else if (new_size < cur_size/2) then
       ! Reduce number of particles
       chance = new_size / real(cur_size, dp)
       do ix = 1, cur_size
          if (pc%rng%uni_01() > chance) call pc%remove_part(ix)
       end do
       call pc%clean_up()
    end if
  end subroutine resize_swarm

  subroutine update_particle_stats(pc, ps, new_ps)
    type(PC_t), intent(inout)         :: pc
    type(part_stats_t), intent(inout) :: ps
    logical, intent(in)               :: new_ps

    integer                           :: ix, num_part
    real(dp)                          :: inv_n_samples
    real(dp)                          :: v(3), v2(3), corr_fac
    real(dp)                          :: coll_rates(PC_max_num_coll)
    type(PC_part_t)                   :: part

    call recenter_swarm(pc)

    if (new_ps) then
       ps%n_samples = 0
       ps%v         = 0.0_dp
       ps%v2        = 0.0_dp
       ps%cov_xv    = 0.0_dp
       ps%coll_rate = 0.0_dp
       ps%i_rate    = 0.0_dp
       ps%a_rate    = 0.0_dp
    end if

    num_part = pc%get_num_sim_part()
    corr_fac = num_part/(num_part-1.0_dp)

    do ix = 1, num_part
       part = pc%get_part(ix)
       ps%n_samples  = ps%n_samples + 1
       inv_n_samples = 1.0_dp / ps%n_samples
       v             = part%v
       v2            = v**2
       ps%v          = ps%v + (v - ps%v) * inv_n_samples
       ! ps%v2_v       = ps%v2_v + (v2 * v - ps%v2_v) * inv_n_samples

       ! This placing is intentional:
       ! http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
       ps%cov_xv     = ps%cov_xv + (v - ps%v) * part%x * corr_fac

       ! ps%cov_v2_xv  = ps%cov_v2_xv + v2 * (v - ps%v) * part%x
       ps%v2         = ps%v2 + (v2 - ps%v2) * inv_n_samples

       call pc%get_coll_rates(sqrt(sum(v2)), coll_rates(1:pc%n_colls))
       ps%coll_rate = ps%coll_rate + (sum(coll_rates(1:pc%n_colls)) - &
            ps%coll_rate) * inv_n_samples
       ps%i_rate = ps%i_rate + (sum(coll_rates(pc%ionization_colls)) - &
            ps%i_rate) * inv_n_samples
       ps%a_rate = ps%a_rate + (sum(coll_rates(pc%attachment_colls)) - &
            ps%a_rate) * inv_n_samples
    end do
  end subroutine update_particle_stats

  subroutine update_td_from_ps(tds, ps)
    use m_units_constants
    type(SWARM_td_t), intent(inout) :: tds(:)
    type(part_stats_t), intent(in)  :: ps

    call update_td(tds(ix_vel_sq), ps%v2)
    call update_td(tds(ix_vel), ps%v)
    call update_td(tds(ix_diffusion), ps%cov_xv / ps%n_samples)
    call update_td(tds(ix_alpha), [ps%i_rate / norm2(ps%v)])
    call update_td(tds(ix_eta), [ps%a_rate / norm2(ps%v)])
    call update_td(tds(ix_coll_rate), [ps%coll_rate])

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

  subroutine new_swarm(pc, swarm_size)
    use m_units_constants
    use m_random
    type(PC_t), intent(inout)       :: pc
    integer, intent(in)             :: swarm_size
    integer                         :: ll
    real(dp)                        :: x(3), v(3), a(3), vel
    real(dp), parameter             :: eV = 1 ! Create 1 eV particles

    call pc%remove_particles()

    do ll = 1, swarm_size
       x   = 0
       a   = SWARM_field%E_vec * UC_elec_q_over_m
       vel = sqrt(2 * eV * UC_elem_charge/UC_elec_mass)
       v   = pc%rng%sphere(vel)
       call pc%create_part(x, v, a, 1.0_dp, 0.0_dp)
    end do
  end subroutine new_swarm

  subroutine SWARM_get_data(pc, swarm_size, tds)
    use iso_fortran_env, only: error_unit
    use m_units_constants

    type(PC_t), intent(inout)       :: pc
    integer, intent(in)             :: swarm_size
    type(SWARM_td_t), intent(inout) :: tds(:)

    integer, parameter  :: n_swarms_min = 10
    integer, parameter  :: n_swarms_max = 10000
    real(dp), parameter :: fac          = 0.5 * UC_elec_mass/UC_elem_charge
    integer             :: n_coll_times
    integer             :: n, n_swarms
    real(dp)            :: tau_coll
    real(dp)            :: growth_rate
    type(part_stats_t)  :: ps

    call create_swarm(pc, tau_coll, swarm_size)
    call update_particle_stats(pc, ps, .true.)

    ! Loop over the swarms until converged
    do n_swarms = 1, n_swarms_max

       ! Estimate the time scale for energy relaxation, given by:
       ! energy / (d/dt energy), in units of the collision time.
       n_coll_times = nint(fac * sum(ps%v2) / &
            abs(dot_product(SWARM_field%E_vec, ps%v) * tau_coll))

       ! For safety, limit n_coll_times to 10 -- 2000
       n_coll_times = max(10, n_coll_times)
       n_coll_times = min(2000, n_coll_times)

       growth_rate  = abs(ps%i_rate - ps%a_rate)

       do n = 1, n_coll_times
          call swarm_advance(pc, tau_coll, swarm_size, growth_rate)
          call update_particle_stats(pc, ps, n == 1)
       end do

       call update_td_from_ps(tds, ps)

       if (n_swarms >= n_swarms_min) then
          ! Check whether the results are accurate enough
          if (check_accuracy(tds)) exit
       end if

       ! Advance over several collisions for the diffusion measurements
       call collapse_swarm(pc)
       call swarm_advance(pc, n_coll_times * tau_coll, &
            swarm_size, growth_rate)
    end do

    if (n_swarms == n_swarms_max + 1) then
       write(error_unit, *) "No convergence in ", n_swarms_max, "iterations"
       error stop
    end if

  end subroutine SWARM_get_data

  !> Check whether the transport data td with estimated deviation dev is
  !> accurate enough according to the requirements in acc
  logical function check_accuracy(tds)
    type(SWARM_td_t), intent(in) :: tds(:) !< The transport data
    real(dp)                     :: stddev, mean, var
    integer                      :: i, n

    check_accuracy = .true.

    do i = 1, SWARM_num_td
       ! Below we use the norm for vector-based transport data, so that
       ! orientation of the field should not matter
       n      = tds(i)%n_measurements
       mean   = norm2(tds(i)%val)             ! Mean
       var    = norm2(tds(i)%var / (n-1)) / n ! Variance of the mean
       stddev = sqrt(var)                     ! Standard deviation of the mean

       if (stddev > tds(i)%rel_acc * mean .and. stddev > tds(i)%abs_acc) then
          check_accuracy = .false.
          exit
       end if
    end do
  end function check_accuracy

  ! Create a swarm that is relaxed to the electric field
  subroutine create_swarm(pc, tau_coll, swarm_size)
    use m_units_constants
    !> Data type which stores particle model
    type(PC_t), intent(inout)       :: pc
    !> Estimate of time between collisions
    real(dp), intent(out)           :: tau_coll
    !> Number of electrons in swarm
    integer, intent(in)             :: swarm_size

    integer, parameter    :: frame_size    = 100
    real(dp), parameter   :: en_eV         = 0.1_dp
    integer, parameter    :: max_its_relax = 500
    integer, parameter    :: min_its_relax = 5
    integer               :: i, ll, cntr
    real(dp)              :: en_hist(frame_size), t_hist(frame_size)
    real(dp)              :: mean_en, correl, stddev, tau
    real(dp), allocatable :: coll_rates(:), tmp_vec(:)


    call new_swarm(pc, swarm_size)

    ! An electron accelerating from zero velocity gains en_eV in this time (in
    ! the absence of a magnetic field). This time step is only used to determine
    ! when the swarm is approximately relaxed to the background field.
    tau = sqrt(0.5_dp * en_eV * UC_elec_mass / UC_elem_charge) / &
         norm2(SWARM_field%E_vec)

    ! Create linear table with unit variance and zero mean
    do i = 1, frame_size
       t_hist(i) = (i - 0.5_dp * (frame_size+1)) * &
            sqrt(12.0_dp / (frame_size + frame_size**2))
    end do

    ! Determine when the mean energy is relaxed, so that we can use this swarm
    do cntr = 1, max_its_relax
       do i = 1, frame_size
          call pc%advance(tau)
          call resize_swarm(pc, swarm_size)
          en_hist(i) = pc%get_mean_energy() / UC_elec_volt
       end do

       mean_en = sum(en_hist) / frame_size
       stddev = sqrt(sum((en_hist - mean_en)**2) / (frame_size-1))

       ! Compute correlation between en_hist and a line
       correl = sum((en_hist - mean_en) * t_hist) / (frame_size * stddev)

       ! If the correlation is sufficiently small, exit
       if (cntr > min_its_relax .and. abs(correl) < 0.25_dp) exit
    end do

    ! Swarm is relaxed, place it at the origin
    call recenter_swarm(pc)

    ! Get current mean collision rates
    allocate(coll_rates(pc%n_colls))
    allocate(tmp_vec(pc%n_colls))
    coll_rates = 0

    do ll = 1, pc%n_part
       call pc%get_coll_rates(norm2(pc%particles(ll)%v), tmp_vec)
       coll_rates = coll_rates + tmp_vec / pc%n_part
    end do

    ! Return estimate of collision time
    tau_coll = 1 / sum(coll_rates)

  end subroutine create_swarm

  subroutine SWARM_print_results(tds)
    use m_units_constants
    type(SWARM_td_t), intent(in) :: tds(:) !< The transport data
    integer                      :: i, i_dim
    real(dp)                     :: fac, tmp, std
    real(dp) :: energy, mu

    i = tds(1)%n_measurements
    fac = sqrt(1.0_dp / (i * (i-1)))

    write(*, "(A24,2E12.4)") "Bz = ", SWARM_field%Bz, 0.0_dp
    write(*, "(A24,2E12.4)") "Ez = ", SWARM_field%Ez, 0.0_dp
    write(*, "(A24,2E12.4)") "Ey = ", SWARM_field%Ey, 0.0_dp
    write(*, "(A24,2E12.4)") "angle = ", SWARM_field%angle_deg, 0.0_dp
    write(*, "(A24,2E12.4)") "omega_c = ", SWARM_field%omega_c, 0.0_dp

    ! mean energy
    tmp    = 0.5_dp * UC_elec_mass / UC_elec_volt
    energy = tmp * sum(tds(ix_vel_sq)%val)
    std    = tmp * sqrt(sum(tds(ix_vel_sq)%var)) * fac
    write(*, "(A24,2E12.4)") "energy = ", energy, std

    ! mobility parallel to E
    mu = -dot_product(tds(ix_vel)%val,  SWARM_field%E_vec) &
         / max(epsilon(1.0_dp), sum(SWARM_field%E_vec**2))
    std = fac * sqrt(dot_product(tds(ix_vel)%var,  SWARM_field%E_vec)) &
         / max(epsilon(1.0_dp), norm2(SWARM_field%E_vec)**1.5_dp)
    write(*, "(A24,2E12.4)") "mu_E = ", mu, std

    ! mobility parallel to B
    if (abs(SWARM_field%Ez) > sqrt(epsilon(1.0_dp))) then
       mu = -tds(ix_vel)%val(3) / SWARM_field%Ez
       std = fac * sqrt(tds(ix_vel)%var(3)) / abs(SWARM_field%Ez)
    else
       mu = 0
       std = 0
    end if
    write(*, "(A24,2E12.4)") "mu_B = ", mu, std

    ! mobility perpendicular to B (y-velocity over Ey)
    if (abs(SWARM_field%Ey) > sqrt(epsilon(1.0_dp))) then
       mu = -tds(ix_vel)%val(2) / SWARM_field%Ey
       std = fac * sqrt(tds(ix_vel)%var(2)) / abs(SWARM_field%Ey)
    else
       mu = 0
       std = 0
    end if
    write(*, "(A24,2E12.4)") "mu_xB = ", mu, std

    ! mobility in ExB-direction (x-velocity over Ey = E_perp)
    if (abs(SWARM_field%Ey) > sqrt(epsilon(1.0_dp))) then
       mu = tds(ix_vel)%val(1) / SWARM_field%Ey
       std = fac * sqrt(tds(ix_vel)%var(1)) / abs(SWARM_field%Ey)
    else
       mu = 0
       std = 0
    end if
    write(*, "(A24,2E12.4)") "mu_ExB = ", mu, std

    ! The other swarm parameters
    do i = 1, SWARM_num_td
       if (tds(i)%n_dim == 1) then
          write(*, "(A24,2E12.4)") trim(tds(i)%description) // " = ", &
               tds(i)%val(1), sqrt(tds(i)%var(1)) * fac
       else

          do i_dim = 1, tds(i)%n_dim
             write(*, "(A20,I0,A,2E12.4)") trim(tds(i)%description) // "_", &
                  i_dim, " = ", tds(i)%val(i_dim), &
                  sqrt(tds(i)%var(i_dim)) * fac
          end do
       end if
    end do

  end subroutine SWARM_print_results

  !> Advance the particle position and velocity over time dt taking into account
  !> a constant electric and magnetic field, using the analytic solution.
  subroutine SWARM_particle_mover_analytic(part, dt)
    use m_units_constants
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
    part%x = part%x + (rc_rot - rc)

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
  subroutine SWARM_particle_mover_simple(part, dt)
    use m_units_constants
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
  subroutine SWARM_particle_mover_boris(part, dt)
    use m_units_constants
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
