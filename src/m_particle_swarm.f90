module m_particle_swarm
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type part_stats_t
     integer            :: n_samples
     real(dp)           :: en
     real(dp)           :: fld
     real(dp)           :: i_rate
     real(dp)           :: a_rate
     real(dp)           :: v(3)
     real(dp)           :: v2
     real(dp)           :: cov_xv(3)
  end type part_stats_t

  integer, parameter :: SWARM_num_td = 6
  integer, parameter :: SWARM_ix_fld = 1
  integer, parameter :: SWARM_ix_en = 2
  integer, parameter :: SWARM_ix_mu = 3
  integer, parameter :: SWARM_ix_D = 4
  integer, parameter :: SWARM_ix_alpha = 5
  integer, parameter :: SWARM_ix_eta = 6

  type SWARM_acc_t
     real(dp) :: energy(2)
     real(dp) :: mobility(2)
     real(dp) :: diffusion(2)
  end type SWARM_acc_t

  character(len=20) :: SWARM_td_names(SWARM_num_td) = &
       [character(len=20) :: "field[V/m]", "en[eV]", &
       "mu[m2/(Vs)]", "Dc[m2/s]", "alpha[1/m]", "eta[1/m]"]

  public :: SWARM_acc_t
  public :: SWARM_num_td
  public :: SWARM_td_names

  ! Public routines from this module
  public :: SWARM_get_data

contains

  subroutine swarm_advance(pc, tau, desired_num_part, growth_rate)
    type(PC_t), intent(inout) :: pc
    real(dp), intent(in)      :: tau
    integer, intent(in)       :: desired_num_part
    integer                   :: n, n_steps
    real(dp)                  :: dt, growth_rate

    ! Sometimes a swarm can rapidly grow or shink in time. Therefore we advance
    ! the swarm in steps, so that we can resize it if necessary.
    n_steps = ceiling(tau * growth_rate / log(2.0_dp))
    n_steps = max(n_steps, 1)
    dt = tau/n_steps

    do n = 1, n_steps
       call pc%advance(dt)
       call resize_swarm(pc, desired_num_part)
    end do
  end subroutine swarm_advance

  subroutine recenter_swarm(pc)
    type(PC_t), intent(inout) :: pc
    real(dp)                  :: sum_x(3), avg_x(3)
    call pc%compute_vector_sum(get_position, sum_x)
    avg_x = sum_x / pc%get_num_sim_part()
    call pc%translate(-avg_x)
  end subroutine recenter_swarm

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

  subroutine get_td_from_ps(ps, fld, td)
    use m_units_constants

    type(part_stats_t), intent(in) :: ps
    real(dp), intent(in)              :: fld
    real(dp), intent(out)             :: td(SWARM_num_td)
    real(dp)                          :: drift_vel

    drift_vel = abs(ps%v(3))
    td(1)     = fld
    td(2)     = 0.5 * UC_elec_mass * ps%v2 / UC_elec_volt ! Energy
    td(3)     = abs(drift_vel / fld)                      ! Mobility
    td(4)     = ps%cov_xv(3) / (ps%n_samples)             ! Long. diffusion
    td(5)     = ps%i_rate / drift_vel                     ! Ionization
    td(6)     = ps%a_rate / drift_vel                     ! Attachment
    ! td(7)     = abs(ps%v2_v(3) / (fld * ps%v2))           ! Energy mobility
    ! td(8)     = ps%cov_v2_xv(1) / (ps%n_samples * ps%v2)  ! Long. energy diffusion
  end subroutine get_td_from_ps

  subroutine new_swarm(pc, swarm_size, fld)
    use m_units_constants
    use m_random
    type(PC_t), intent(inout) :: pc
    integer, intent(in)       :: swarm_size
    real(dp), intent(in)      :: fld
    integer                   :: ll
    real(dp)                  :: x(3), v(3), a(3), vel
    real(dp), parameter       :: eV = 1 ! Create 1 eV particles

    call pc%remove_particles()

    do ll = 1, swarm_size
       x = 0
       a = (/0.0_dp, 0.0_dp, fld/) * UC_elec_charge / UC_elec_mass
       vel = sqrt(2 * eV * UC_elem_charge/UC_elec_mass)
       v = pc%rng%sphere(vel)
       call pc%create_part(x, v, a, 1.0_dp, 0.0_dp)
    end do
  end subroutine new_swarm

  subroutine SWARM_get_data(pc, fld, swarm_size, n_swarms_min, &
       acc, td, td_dev)
    use m_units_constants

    type(PC_t), intent(inout)     :: pc
    real(dp), intent(in)          :: fld
    integer, intent(in)           :: swarm_size, n_swarms_min
    type(SWARM_acc_t), intent(in) :: acc
    real(dp), intent(inout)       :: td(SWARM_num_td)
    real(dp), intent(out)         :: td_dev(SWARM_num_td)

    real(dp), parameter :: fac = 0.5 * UC_elec_mass/UC_elem_charge
    integer             :: n_coll_times
    integer             :: n, n_swarms
    real(dp)            :: tau_coll, weight
    real(dp)            :: td_prev(SWARM_num_td)
    real(dp)            :: abs_acc(SWARM_num_td)
    real(dp)            :: rel_acc(SWARM_num_td)
    real(dp)            :: growth_rate
    logical             :: accurate
    type(part_stats_t)  :: ps

    n_swarms = 0
    td_dev   = 0.0_dp
    td_prev  = 0.0_dp
    accurate = .false.

    ! Set accuracy requirements
    abs_acc = huge(1.0_dp)
    rel_acc = 0
    rel_acc(SWARM_ix_en) = acc%energy(1)
    abs_acc(SWARM_ix_en) = acc%energy(2)
    rel_acc(SWARM_ix_mu) = acc%mobility(1)
    abs_acc(SWARM_ix_mu) = acc%mobility(2)
    rel_acc(SWARM_ix_D)  = acc%diffusion(1)
    abs_acc(SWARM_ix_D)  = acc%diffusion(2)

    call create_swarm(pc, fld, tau_coll, swarm_size)

    ps%fld = fld
    call update_particle_stats(pc, ps, .true.)

    ! Loop over the swarms until converged
    do while (.not. accurate .or. n_swarms < n_swarms_min)

       ! Estimate the time scale for energy relaxation, given by:
       ! energy / (d/dt energy), in units of the collision time.
       n_coll_times = nint(fac * ps%v2 / abs(ps%fld * ps%v(3) * tau_coll))

       ! For safety, limit n_coll_times to 10 -- 2000
       n_coll_times = max(10, n_coll_times)
       n_coll_times = min(2000, n_coll_times)

       growth_rate = abs(ps%i_rate - ps%a_rate)
       do n = 1, n_coll_times
          call swarm_advance(pc, tau_coll, swarm_size, growth_rate)
          call update_particle_stats(pc, ps, .false.)
       end do

       n_swarms = n_swarms + 1
       call get_td_from_ps(ps, fld, td)
       if (n_swarms == 1) then
          td_prev  = td
       else
          weight = 1.0_dp / min(n_swarms-1, n_swarms_min)
          td_dev   = (1 - weight) * td_dev + &
               sqrt(real(n_swarms, dp)) * abs(td - td_prev) * weight
          td_prev  = td
          accurate = all(td_dev < abs_acc .or. td_dev < rel_acc * td)
       end if

       ! Advance over several collisions for the diffusion measurements
       call collapse_swarm(pc)
       call swarm_advance(pc, n_coll_times * tau_coll, &
            swarm_size, growth_rate)
    end do
  end subroutine SWARM_get_data

  ! Create a swarm that is relaxed to the electric field
  subroutine create_swarm(pc, fld, tau_coll, swarm_size)
    use m_units_constants
    !> Data type which stores particle model
    type(PC_t), intent(inout) :: pc
    !> Electric field in which the swarm will propagate
    real(dp), intent(in)      :: fld
    !> Estimate of time between collisions
    real(dp), intent(out)     :: tau_coll
    !> Number of electrons in swarm
    integer, intent(in)       :: swarm_size

    integer, parameter    :: frame_size    = 100
    real(dp), parameter   :: en_eV         = 0.1_dp
    integer, parameter    :: max_its_relax = 500
    integer, parameter    :: min_its_relax = 5
    integer               :: i, ll, cntr
    real(dp)              :: en_hist(frame_size), t_hist(frame_size)
    real(dp)              :: mean_en, correl, stddev, tau
    real(dp), allocatable :: coll_rates(:), tmp_vec(:)


    call new_swarm(pc, swarm_size, fld)

    ! An electron accelerating from zero velocity gains en_eV in this time. This
    ! time step is only used to determine when the swarm is approximately
    ! relaxed to the background field.
    tau = sqrt(0.5_dp * en_eV * UC_elec_mass / UC_elem_charge) / abs(fld)

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

end module m_particle_swarm
