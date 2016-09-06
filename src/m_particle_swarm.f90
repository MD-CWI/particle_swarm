module m_particle_swarm
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  !> Index of energy in transport data list
  integer, parameter :: i_energy     = 1
  !> Index of mobility(3) in transport data list
  integer, parameter :: i_mobility   = i_energy + 1
  !> Index of diffusion(2) in transport data list
  integer, parameter :: i_diffusion  = i_mobility + 3
  !> Index of ionization coeff. in transport data list
  integer, parameter :: i_ionization = i_diffusion + 2
  !> Index of attachment coeff. in transport data list
  integer, parameter :: i_attachment = i_ionization + 1

  !> The number of transport data fields
  integer, parameter :: SWARM_num_td = i_attachment

  character(len=20), parameter :: SWARM_td_names(SWARM_num_td) = &
       [character(len=20) :: "energy (eV)", "mu_x", "mu_y", "mu_z", &
       "D_L", "D_T", "alpha", "eta"]

  !> Type storing accuracy requirements for a swarm
  type SWARM_acc_t
     real(dp) :: relative(SWARM_num_td) = 0.0_dp
     real(dp) :: absolute(SWARM_num_td) = huge(1.0_dp)
  end type SWARM_acc_t

  !> Type storing the field configuration for a swarm
  type SWARM_field_t
     real(dp) :: B_vec(3) !< Magnetic field vector
     real(dp) :: E_vec(3) !< Electric field vector
     real(dp) :: Bz       !< Magnetic field along z-axis (T)
     real(dp) :: Ez       !< Electric field along z-axis (V/m)
     real(dp) :: Ey       !< Electric field along y-axis (V/m)
  end type SWARM_field_t

  !> Type to collect particle statistics / properties
  type part_stats_t
     type(SWARM_field_t) :: field     !< The field configuration
     integer             :: n_samples !< Number of samples
     real(dp)            :: i_rate    !< Ionization rate
     real(dp)            :: a_rate    !< Attachment rate
     real(dp)            :: v(3)      !< Mean velocity
     real(dp)            :: v2        !< Mean square velocity
     real(dp)            :: cov_xv(3) !< Covariance x,v
  end type part_stats_t

  type(SWARM_field_t), protected :: SWARM_field
  real(dp), parameter :: SWARM_omega_unitvec(3) = [0.0_dp, 0.0_dp, 1.0_dp]
  real(dp) :: SWARM_omega_c
  real(dp) :: SWARM_plasma_drift_velocity(3)

  public :: i_energy
  public :: i_mobility
  public :: i_diffusion
  public :: i_ionization
  public :: i_attachment
  public :: SWARM_num_td
  public :: SWARM_acc_t
  public :: SWARM_field_t

  ! Public routines from this module
  public :: SWARM_get_data
  public :: SWARM_print_results

contains

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
    real(dp)                          :: v(3), v2, corr_fac
    real(dp)                          :: coll_rates(PC_max_num_coll)
    type(PC_part_t)                   :: part

    call recenter_swarm(pc)

    if (new_ps) then
       ps%n_samples = 0
       ps%v         = 0.0_dp
       ps%v2        = 0.0_dp
       ps%cov_xv    = 0.0_dp
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
       v2            = sum(v**2)
       ps%v          = ps%v + (v - ps%v) * inv_n_samples
       ! ps%v2_v       = ps%v2_v + (v2 * v - ps%v2_v) * inv_n_samples

       ! This placing is intentional:
       ! http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
       ps%cov_xv     = ps%cov_xv + (v - ps%v) * part%x * corr_fac

       ! ps%cov_v2_xv  = ps%cov_v2_xv + v2 * (v - ps%v) * part%x
       ps%v2         = ps%v2 + (v2 - ps%v2) * inv_n_samples
       call pc%get_coll_rates(sqrt(v2), coll_rates(1:pc%n_colls))
       ps%i_rate = ps%i_rate + (sum(coll_rates(pc%ionization_colls)) - &
            ps%i_rate) * inv_n_samples
       ps%a_rate = ps%a_rate + (sum(coll_rates(pc%attachment_colls)) - &
            ps%a_rate) * inv_n_samples
    end do
  end subroutine update_particle_stats

  subroutine get_td_from_ps(ps, td)
    use m_units_constants

    type(part_stats_t), intent(in)     :: ps
    real(dp), intent(out) :: td(SWARM_num_td)

    td(i_energy)                = 0.5 * UC_elec_mass * ps%v2 / UC_elec_volt
    td(i_mobility:i_mobility+2) = ps%v / norm2(SWARM_field%E_vec)             ! Mobility
    td(i_diffusion)             = ps%cov_xv(3) / ps%n_samples                 ! Long. diffusion
    td(i_diffusion+1)           = 0.5_dp * sum(ps%cov_xv(1:2)) / ps%n_samples ! Trans. diffusion
    td(i_ionization)            = ps%i_rate / norm2(ps%v)
    td(i_attachment)            = ps%a_rate / norm2(ps%v)
    ! td(7)     = abs(ps%v2_v(3) / (fld * ps%v2))           ! Energy mobility
    ! td(8)     = ps%cov_v2_xv(1) / (ps%n_samples * ps%v2)  ! Long. energy diffusion
  end subroutine get_td_from_ps

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
       a = SWARM_field%E_vec * UC_elec_q_over_m
       vel = sqrt(2 * eV * UC_elem_charge/UC_elec_mass)
       v   = pc%rng%sphere(vel)
       call pc%create_part(x, v, a, 1.0_dp, 0.0_dp)
    end do
  end subroutine new_swarm

  subroutine SWARM_get_data(pc, field, swarm_size, acc, td, td_dev)
    use m_units_constants

    type(PC_t), intent(inout)       :: pc
    type(SWARM_field_t), intent(in) :: field
    integer, intent(in)             :: swarm_size
    type(SWARM_acc_t), intent(in)   :: acc
    real(dp), intent(out)           :: td(SWARM_num_td)
    real(dp), intent(out)           :: td_dev(SWARM_num_td)

    integer, parameter  :: n_swarms_min = 10
    real(dp), parameter :: fac          = 0.5 * UC_elec_mass/UC_elem_charge
    integer             :: n_coll_times
    integer             :: n, n_swarms
    real(dp)            :: tau_coll
    real(dp)            :: td_prev(SWARM_num_td)
    real(dp)            :: growth_rate
    integer             :: n_accurate
    type(part_stats_t)  :: ps

    n_swarms = 0
    td_prev  = 0.0_dp
    n_accurate = 0

    SWARM_field = field

    if (abs(field%Bz) > 0) then
       ! Check whether the magnetic field is large enough
       SWARM_omega_c = abs(UC_elem_charge * field%Bz / UC_elec_mass)
       SWARM_plasma_drift_velocity = &
            [field%Ey * field%Bz, 0.0_dp, 0.0_dp] / field%Bz**2

       pc%particle_mover => advance_particle_analytic
    end if

    call create_swarm(pc, tau_coll, swarm_size)

    call update_particle_stats(pc, ps, .true.)

    ! Loop over the swarms until converged
    do while (n_swarms < n_swarms_min .or. n_accurate < n_swarms_min)

       ! Estimate the time scale for energy relaxation, given by:
       ! energy / (d/dt energy), in units of the collision time.
       n_coll_times = nint(fac * ps%v2 / abs(dot_product(SWARM_field%E_vec, ps%v) * tau_coll))

       ! For safety, limit n_coll_times to 10 -- 2000
       n_coll_times = max(10, n_coll_times)
       n_coll_times = min(2000, n_coll_times)

       growth_rate = abs(ps%i_rate - ps%a_rate)
       do n = 1, n_coll_times
          call swarm_advance(pc, tau_coll, swarm_size, growth_rate)
          call update_particle_stats(pc, ps, .false.)
       end do

       n_swarms = n_swarms + 1
       call get_td_from_ps(ps, td)

       if (n_swarms == 1) then
          td_prev = td
       else
          ! Using the current and previous transport data, estimate the
          ! standard deviation (or the error)
          td_dev = sqrt(real(n_swarms, dp)) * abs(td - td_prev)

          ! Check whether the results are accurate enough
          if (check_accuracy(td, td_dev, acc)) n_accurate = n_accurate + 1

          td_prev = td
       end if

       ! Advance over several collisions for the diffusion measurements
       call collapse_swarm(pc)
       call swarm_advance(pc, n_coll_times * tau_coll, &
            swarm_size, growth_rate)
    end do
  end subroutine SWARM_get_data

  !> Check whether the transport data td with estimated deviation dev is
  !> accurate enough according to the requirements in acc
  logical function check_accuracy(td, dev, acc)
    real(dp), intent(in)          :: td(SWARM_num_td)  !< The transport data
    real(dp), intent(in)          :: dev(SWARM_num_td) !< The estimated error in td
    type(SWARM_acc_t), intent(in) :: acc               !< The accuracy requirements
    integer                       :: i

    check_accuracy = .true.

    do i = 1, SWARM_num_td
       if (dev(i) > acc%relative(i) * abs(td(i)) .and. &
            dev(i) > acc%absolute(i)) then
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

    print *, tau

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

  subroutine SWARM_print_results(td, dev)
    real(dp), intent(in)            :: td(SWARM_num_td)  !< The transport data
    real(dp), intent(in)            :: dev(SWARM_num_td) !< The estimated error in td
    integer                         :: i

    write(*, "(A20,E10.3)") " Bz (T)             ", SWARM_field%Bz
    write(*, "(A20,E10.3)") " Ez (V/m)           ", SWARM_field%Ez
    write(*, "(A20,E10.3)") " Ey (V/m)           ", SWARM_field%Ey

    do i = 1, SWARM_num_td
       write(*, "(A20,E10.3,E9.2)") " " // SWARM_td_names(i), td(i), dev(i)
    end do

  end subroutine SWARM_print_results

  !> Advance the particle position and velocity over time dt taking into account
  !> a constant electric and magnetic field
  subroutine advance_particle_analytic(part, dt)
    use m_units_constants
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: dt

    real(dp) :: rc(3), rc_rot(3), theta

    ! First the motion perpendicular to the magnetic field

    ! Subtract the plasma drift velocity
    part%v = part%v - SWARM_plasma_drift_velocity

    ! Determine rc, the vector pointing from the center of gyration
    rc = cross_product(part%v, SWARM_omega_unitvec) / SWARM_omega_c

    ! Determine rotation angle
    theta = dt * SWARM_omega_c

    ! Rotate position and velocity
    part%v = rotate_around_axis(part%v, SWARM_omega_unitvec, theta)
    rc_rot = rotate_around_axis(rc, SWARM_omega_unitvec, theta)

    ! Update the position with the change in guiding center
    part%x = part%x + (rc_rot - rc)

    ! Add back the plasma drift velocity
    part%v = part%v + SWARM_plasma_drift_velocity

    ! Now the motion parallel to the magnetic field
    part%v(3) = part%v(3) + 0.5_dp * dt * UC_elec_q_over_m * SWARM_field%Ez
    part%x(3) = part%x(3) + dt * part%v(3)
    part%v(3) = part%v(3) + 0.5_dp * dt * UC_elec_q_over_m * SWARM_field%Ez

    ! Update time left
    part%t_left = part%t_left - dt
  end subroutine advance_particle_analytic

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
