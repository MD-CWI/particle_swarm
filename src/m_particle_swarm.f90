module m_particle_swarm
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type PM_part_stats_t
     integer            :: n_samples
     real(dp)           :: en
     real(dp)           :: fld
     real(dp)           :: i_rate
     real(dp)           :: a_rate
     real(dp)           :: v(3)
     real(dp)           :: v2
     real(dp)           :: v2_v(3)
     real(dp)           :: cov_xv(3)
     real(dp)           :: cov_v2_xv(3)
  end type PM_part_stats_t

  integer, parameter :: PM_num_td = 8
  character(len=20) :: PM_td_names(PM_num_td) = (/ &
       "field[V/m]    ", &
       "en[eV]        ", &
       "mu[m2/(Vs)]   ", &
       "Dc[m2/s]      ", &
       "alpha[1/m]    ", &
       "eta[1/m]      ", &
       "mu_en[m2/(Vs)]", &
       "Dc_en[m2/s]   " /)


  integer :: PM_num_coll = 0
  integer, allocatable :: PM_ionization_colls(:)
  integer, allocatable :: PM_attachment_colls(:)

  public :: PM_num_td
  public :: PM_td_names
  public :: PM_part_stats_t

  ! Public routines from particle core module
  public :: PC_get_num_sim_part
  public :: PC_advance
  public :: PC_clean_up_dead
  public :: PC_get_max_coll_rate
  public :: PC_reset_coll_count

  ! Public routines from this module
  public :: PM_advance
  public :: PM_get_avg_growth_rate
  public :: PM_initialize
  public :: PM_update_particle_stats
  public :: PM_resize_swarm
  public :: PM_recenter_swarm
  public :: PM_new_swarm
  public :: PM_get_mean_energy
  public :: PM_get_td_from_ps
  public :: PM_collapse_swarm
  public :: PM_get_avg_coll_rate

contains

  subroutine PM_advance(tau, desired_num_part, avg_growth_rate)
    real(dp), intent(in) :: tau
    integer, intent(in) :: desired_num_part
    integer :: n, n_steps
    real(dp) :: dt, avg_growth_rate

    n_steps = ceiling(log(1.5_dp) * tau * avg_growth_rate)
    n_steps = max(n_steps, 1)
    dt = tau/n_steps

    do n = 1, n_steps
       call PC_advance(dt)
       call PM_resize_swarm(desired_num_part)
    end do
  end subroutine PM_advance

  subroutine recenter_swarm(pc)
    type(PC_t), intent(inout) :: pc
    real(dp) :: sum_x(3), avg_x(3)
    call pc%compute_vector_sum(get_position, sum_x)
    avg_x = sum_x / pc%get_num_sim_part()
    call pc%translate(-avg_x)
  end subroutine recenter_swarm

  subroutine PM_collapse_swarm(pc)
    type(PC_t), intent(in) :: pc
    call pc%loop_iopart(reset_part_x)
  end subroutine PM_collapse_swarm

  subroutine get_position(part, x)
    type(PC_part_t), intent(in) :: part
    real(dp), intent(out)       :: x(:)
    x = part%x
  end subroutine get_position

  subroutine reset_part_x(part)
    type(PC_part_t), intent(inout) :: part
    part%x = 0.0_dp
  end subroutine reset_part_x

  subroutine PM_resize_swarm(new_size)
    use m_random
    integer, intent(in) :: new_size

    integer :: ix, cur_size
    real(dp) :: chance, rand_val
    type(PC_part_t) :: part

    call PC_clean_up_dead()
    cur_size = PC_get_num_sim_part()

    if (new_size >= 2 * cur_size) then
       do
          ! Double particles
          do ix = 1, cur_size
             call PC_get_part(ix, part)
             call PC_add_part(part)
          end do
          cur_size = cur_size * 2
          if (new_size < 2 * cur_size) exit
       end do
    else if (new_size < cur_size/2) then
       ! Reduce number of particles
       part%live = .false.
       chance = new_size / real(cur_size, dp)
       do ix = 1, cur_size
          call random_number(rand_val)
          if (rand_val > chance) call PC_set_part(ix, part)
       end do
       call PC_clean_up_dead()
    end if
  end subroutine PM_resize_swarm

  subroutine get_ionization_attachment_rate(velocity, two_rates)
    real(dp), intent(in)          :: velocity
    real(dp), intent(out) :: two_rates(2)
    real(dp)                      :: coll_rates(PM_num_coll)
    call PC_get_coll_rates(velocity, coll_rates)
    two_rates(1) = sum(coll_rates(PM_ionization_colls))
    two_rates(2) = sum(coll_rates(PM_attachment_colls))
  end subroutine get_ionization_attachment_rate

  real(dp) function PM_get_avg_coll_rate()
    real(dp) :: sum_rates(PM_num_coll)
    call PC_compute_vector_sum(get_coll_rates, sum_rates)
    PM_get_avg_coll_rate = sum(sum_rates) / PC_get_num_sim_part()
  end function PM_get_avg_coll_rate

  real(dp) function PM_get_avg_growth_rate()
    real(dp) :: sum_rates(PM_num_coll)
    call PC_compute_vector_sum(get_coll_rates, sum_rates)
    PM_get_avg_growth_rate = (sum(sum_rates(PM_ionization_colls)) &
         - sum(sum_rates(PM_attachment_colls))) / PC_get_num_sim_part()
  end function PM_get_avg_growth_rate

  subroutine get_coll_rates(part, rates)
    type(PC_part_t), intent(in) :: part
    real(dp), intent(out) :: rates(:)
    call PC_get_coll_rates(norm2(part%v), rates)
  end subroutine get_coll_rates

  subroutine PM_update_particle_stats(ps, new_ps)
    type(PM_part_stats_t), intent(inout) :: ps
    logical, intent(in)                  :: new_ps

    integer                              :: ix, num_part
    real(dp)                             :: inv_n_samples
    real(dp) :: v(3), v2, two_rates(2)
    type(PC_part_t)                      :: part

    ! Recenter particles
    call PM_recenter_swarm()

    if (new_ps) then
       ps%n_samples = 0
       ps%v         = 0.0_dp
       ps%v2        = 0.0_dp
       ps%v2_v      = 0.0_dp
       ps%cov_xv    = 0.0_dp
       ps%cov_v2_xv = 0.0_dp
       ps%i_rate    = 0.0_dp
       ps%a_rate    = 0.0_dp
    end if

    num_part = PC_get_num_sim_part()

    do ix = 1, num_part
       call PC_get_part(ix, part)
       ps%n_samples  = ps%n_samples + 1
       inv_n_samples = 1.0_dp / ps%n_samples
       v             = part%v
       v2            = sum(v**2)
       ps%v          = ps%v + (v - ps%v) * inv_n_samples
       ps%v2_v       = ps%v2_v + (v2 * v - ps%v2_v) * inv_n_samples

       ! This placing is intentional:
       ! http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
       ps%cov_xv     = ps%cov_xv + (v - ps%v) * part%x
       ps%cov_v2_xv  = ps%cov_v2_xv + v2 * (v - ps%v) * part%x
       ps%v2         = ps%v2 + (v2 - ps%v2) * inv_n_samples

       call get_ionization_attachment_rate(norm2(v), two_rates)
       ps%i_rate     = ps%i_rate + (two_rates(1) - ps%i_rate) * inv_n_samples
       ps%a_rate     = ps%a_rate + (two_rates(2) - ps%a_rate) * inv_n_samples
    end do

  end subroutine PM_update_particle_stats

  real(dp) function PM_get_mean_energy()
    real(dp)        :: sum_en
    call PC_compute_sum(get_energy, sum_en)
    PM_get_mean_energy = sum_en / PC_get_num_sim_part()
  end function PM_get_mean_energy

  subroutine get_energy(part, energy)
    type(PC_part_t), intent(in) :: part
    real(dp), intent(out)       :: energy
    energy = PC_vel_to_en(norm2(part%v))
  end subroutine get_energy

  subroutine PM_get_td_from_ps(ps, fld, td)
    use m_units_constants

    type(PM_part_stats_t), intent(in) :: ps
    real(dp), intent(in)              :: fld
    real(dp), intent(out)             :: td(PM_num_td)
    real(dp)                          :: v_norm

    v_norm = norm2(ps%v)
    td(1) = fld
    td(2) = 0.5 * UC_elec_mass * ps%v2 / UC_elec_volt ! Energy
    td(3) = abs(ps%v(3) / fld)                        ! Mobility
    td(4) = ps%cov_xv(1) / (ps%n_samples)             ! Long. diffusion
    td(5) = ps%i_rate / v_norm                        ! Ionization
    td(6) = ps%a_rate / v_norm                        ! Attachment
    td(7) = abs(ps%v2_v(3) / (fld * ps%v2))           ! Energy mobility
    td(8) = ps%cov_v2_xv(1) / (ps%n_samples * ps%v2)  ! Long. energy diffusion
  end subroutine PM_get_td_from_ps

  subroutine new_swarm(pc, swarm_size, fld)
    use m_units_constants
    use m_random

    integer, intent(in)  :: swarm_size
    real(dp), intent(in) :: fld, init_eV
    integer              :: ll
    real(dp)             :: sigma, mass, x(3), v(3), a(3)

    call pc%reset()
    ! mass = pc%get_mass()
    ! sigma = init_eV * UC_elec_volt / 3
    ! sigma = PC_en_to_vel(sigma, mass)

    do ll = 1, swarmSize
       x = 0
       v = 0
       a = (/0.0_dp, 0.0_dp, fld/) * UC_elec_charge / UC_elec_mass
       call pc%create_part(x, v, a, 1.0_dp, 0.0_dp)
    end do

    call recenter_swarm(pc)
  end subroutine New_swarm

  subroutine PM_get_swarm_data(pmodel, fld, swarm_data)
    type(PC_t), intent(inout) :: pc
    real(dp), intent(in) :: fld
    real(dp), intent(inout) :: swarm_data(:)

    integer :: n_swarm
    real(dp)          :: td_dev(PM_num_td), td_prev(PM_num_td)
    real(dp)          :: td(PM_num_td)
    real(dp)          :: abs_acc(PM_num_td), rel_acc(PM_num_td)

    n_swarms = 0
    td_dev   = 0.0_dp
    td_prev  = 0.0_dp
    accurate = .false.
    print *, "Simulating for fld:", fld, "V/m"

    ! Determine a 'good' initial swarm size and the time until energy
    ! relaxation.
    call create_swarm(pc, init_eV, fld, tau_swarm, tau_coll, growth_rate, swarm_size)
    print *, "Swarm relaxation time", tau_swarm
    print *, "Collision time", tau_coll
    call PM_update_particle_stats(ps, .true.)

    ! Loop over the swarms until converged
    do while (.not. accurate .or. n_swarms < n_swarms_min)
       do mm = 1, max(1, nint(tau_swarm / tau_coll))
          call PM_advance(tau_coll, swarm_size, growth_rate)                  ! Advance the swarm in time
          call PC_clean_up_dead()
          call PM_update_particle_stats(ps, .false.)
       end do

       n_swarms = n_swarms + 1
       call PM_get_td_from_ps(ps, fld, td)
       if (n_swarms == 1) then
          td_prev  = td
       else
          weight = 1.0_dp / min(n_swarms-1, n_swarms_min)
          td_dev   = (1 - weight) * td_dev + &
               sqrt(real(n_swarms, dp)) * abs(td - td_prev) * weight
          td_prev  = td
          accurate = all(td_dev < abs_acc .or. td_dev < rel_acc * td)
       end if

       write(*, fmt="(A)", advance="no") "."
       ! print *, td_dev < abs_acc .or. td_dev < rel_acc * td
       call PM_collapse_swarm()

       ! Advance over 10 collisions for the diffusion measurements
       call PM_advance(10 * tau_coll, swarm_size, growth_rate)
    end do
  end subroutine PM_get_swarm_data

  ! Create a swarm that is relaxed to the electric field
  subroutine create_swarm(energyEv, fld, tau_swarm, tau_coll, growth_rate, swarm_size)
    use m_units_constants
    real(dp), intent(in)  :: fld, energyEv
    real(dp), intent(out) :: tau_swarm, tau_coll, growth_rate
    integer, intent(in)   :: swarm_size

    integer, parameter    :: frame_size = 10
    integer :: i
    real(dp)              :: dev, nu_hist(frame_size)
    real(dp) :: mean_nu, stddev, tau

    call PM_new_swarm(swarm_size, fld, energyEv)

    ! Timestep: an electron accelerating from zero velocity gains 1 eV in this time.
    tau      = 1.0_dp / (sqrt(0.5_dp * UC_elem_charge / UC_elec_mass) * abs(fld))
    tau_swarm = 0.0_dp

    ! Determine when collision rate is relaxed
    do
       do i = 1, frame_size
          call PM_advance(tau, swarm_size, 0.0_dp)
          nu_hist(i) = PM_get_avg_coll_rate()
       end do
       tau_swarm = tau_swarm + tau * frame_size

       mean_nu = sum(nu_hist)/frame_size
       stddev = sqrt(sum((nu_hist(2:) - mean_nu)**2) / frame_size)
       ! Now see whether nu_hist is changing more than stddev in time
       dev    = abs(nu_hist(1) - nu_hist(frame_size)) / stddev
       print *, "Relaxation:", dev

       if (dev < 1.0_dp) exit
    end do

    ! call PM_recenter_swarm()
    tau_coll = 1 / mean_nu
    ! Estimate growth rate
    growth_rate = PM_get_avg_growth_rate()

  end subroutine create_swarm


end module m_particle_swarm
