module m_particle_swarm
  use m_particle_core
  use iso_fortran_env, only: int64

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  !> The axis of gyration for electrons points in the z-direction (just as the
  !> magnetic field)
  real(dp), parameter :: SWARM_omega_unitvec(3) = [0.0_dp, 0.0_dp, 1.0_dp]

  ! These are the basic swarm parameters that we measure, but some others that
  ! can be derived from these (e.g., the mean energy or mobility) are printed as
  ! well.
  integer, parameter :: SWARM_num_td = 10
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
  integer, allocatable :: ix_rates(:)

  !> Indices of ionization collisions
  integer, allocatable :: ionization_colls(:)

  !> Indices of attachment collisions
  integer, allocatable :: attachment_colls(:)

  !> Number of measurements per collision time
  real(dp), parameter :: measurements_per_collision = 1.0_dp

  !> Type for storing transport data
  type SWARM_td_t
     character(len=40)     :: description    !< Description
     character(len=40)     :: unit           !< Unit
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
     integer(int64)      :: n_samples   !< Number of samples
     real(dp)            :: coll_rate   !< Collision rate
     real(dp)            :: i_rate      !< Ionization rate
     real(dp)            :: a_rate      !< Attachment rate
     real(dp), allocatable :: rates(:)  !< All collision rates
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
    type(PC_t), intent(in)                       :: pc
    type(CFG_t), intent(inout)                   :: cfg
    type(SWARM_td_t), allocatable, intent(inout) :: tds(:)
    type(SWARM_field_t), intent(in)              :: field
    real(dp)                                     :: rel_abs_acc(2)
    integer                                      :: n, i_i, i_a
    character(len=40)                            :: name, energy_text, gas

    SWARM_field = field

    allocate(tds(SWARM_num_td + pc%n_colls))
    call init_td(tds(ix_flux_v2), 3, "Velocity squared", "m2/s2")
    call init_td(tds(ix_flux_v), 3, "Drift velocity", "m/s")
    call init_td(tds(ix_bulk_v), 3, "Bulk drift velocity", "m/s")
    call init_td(tds(ix_flux_diff), 3, "Diffusion coefficient", "m2/s")
    call init_td(tds(ix_bulk_diff), 3, "Bulk diffusion coef.", "m2/s")
    call init_td(tds(ix_alpha), 1, "Townsend ioniz. coef. alpha", "1/m")
    call init_td(tds(ix_eta), 1, "Townsend attach. coef. eta", "1/m")
    call init_td(tds(ix_coll_rate), 1, "Total collision rate", "1/s")
    call init_td(tds(ix_ionization), 1, "Total ionization freq.", "1/s")
    call init_td(tds(ix_attachment), 1, "Total attachment freq.", "1/s")

    allocate(ix_rates(pc%n_colls))
    do n = 1, pc%n_colls
       ix_rates(n) = SWARM_num_td + n

       ! Get description of energy threshold
       write(energy_text, '(F10.2,A)') pc%cross_secs(n)%min_energy, " eV"
       energy_text = adjustl(energy_text)
       gas = trim(pc%cross_secs(n)%gas_name)

       select case (pc%cross_secs(n)%coll%type)
       case (CS_ionize_t)
          write(name, '(A,I0,A,A)') "C", n, " " // trim(gas) // &
               " Ionization ", trim(energy_text)
       case (CS_attach_t)
          write(name, '(A,I0,A,A)') "C", n, " " // trim(gas) // &
               " Attachment"
       case (CS_excite_t)
          write(name, '(A,I0,A,A,A)') "C", n, " " // trim(gas) // &
               " Excitation ", trim(energy_text)
       case (CS_elastic_t)
          write(name, '(A,I0,A,A)') "C", n, " " // trim(gas) // &
               " Elastic"
       end select
       call init_td(tds(ix_rates(n)), 1, name, "1/s")
    end do

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
       else if (pc%colls(n)%type == CS_attach_t) then
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

  subroutine resize_swarm(pc, goal_size)
    use m_random
    use m_units_constants
    type(PC_t), intent(inout) :: pc
    !> Desired size of the swarm
    integer, intent(in)       :: goal_size
    integer                   :: ix, ix_group, cur_size, group_size
    real(dp)                  :: chance
    type(PC_part_t)           :: part

    cur_size = pc%get_num_sim_part()

    if (cur_size <= 1 .or. goal_size <= 1) &
         error stop "both goal_size and current size should be > 1"

    if (cur_size <= goal_size/sqrt(2.0_dp)) then
       ! When the swarm decreases in size, aim for a size between
       ! goal_size/sqrt(2) and sqrt(2) * goal_size
       do
          ! Double particles
          do ix = 1, cur_size
             part = pc%particles(ix)
             call pc%add_part(part)
          end do
          cur_size = cur_size * 2
          if (cur_size >= goal_size) exit
       end do
    else if (cur_size >= goal_size * 2) then
       ! Remove particles in stratified way, by first putting them in groups,
       ! from which we remove one particle. Keep size (roughly) between
       ! goal_size and 2 * goal_size.
       chance = goal_size / real(cur_size, dp)
       group_size = max(2, nint(1/chance))

       ! Sort by energy
       call pc%sort(get_v2)

       ! Loop over groups
       do ix_group = 1, cur_size/group_size
          ! First index in group
          ix = (ix_group - 1) * group_size + 1
          ! Add number between 0 and group_size - 1
          ix = ix + floor(pc%rng%unif_01() * group_size)
          call pc%remove_part(ix)
       end do

       ! Handle particles outside last group
       do ix = (cur_size/group_size)*group_size+1, cur_size
          if (pc%rng%unif_01() > 1.0_dp/group_size) call pc%remove_part(ix)
       end do
       call pc%clean_up()
    end if
  end subroutine resize_swarm

  pure real(dp) function get_v2(part)
    type(PC_part_t), intent(in) :: part
    get_v2 = sum(part%v**2)
  end function get_v2

  subroutine initialize_particle_stats(ps, pc)
    type(part_stats_t), intent(inout) :: ps
    type(pc_t), intent(inout)         :: pc
    integer                           :: i, num_part

    if (.not. allocated(ps%rates)) then
       allocate(ps%rates(pc%n_colls))
    else if (size(ps%rates) /= pc%n_colls) then
       deallocate(ps%rates)
       allocate(ps%rates(pc%n_colls))
    end if

    ps%n_samples = 0
    ps%flux_v    = 0.0_dp
    ps%bulk_v    = 0.0_dp
    ps%bulk_dif  = 0.0_dp
    ps%flux_v2   = 0.0_dp
    ps%cov_xv    = 0.0_dp
    ps%coll_rate = 0.0_dp
    ps%i_rate    = 0.0_dp
    ps%a_rate    = 0.0_dp
    ps%rates(:)  = 0.0_dp

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
    real(dp)                          :: fac, inv_n_samples
    real(dp)                          :: x(3), v(3), corr_fac, sum_x2(3)
    real(dp)                          :: sum_cov(3), sum_v2(3)
    real(dp)                          :: sum_v(3), sum_x(3)
    real(dp)                          :: mean_v2(3), mean_x2(3), mean_v(3)
    real(dp)                          :: mean_x(3), bulk_v(3), bulk_dif(3)
    real(dp)                          :: rates(pc%n_colls)
    real(dp)                          :: sum_rates(pc%n_colls)

    num_part  = pc%get_num_sim_part()
    corr_fac  = num_part/(num_part-1.0_dp)
    sum_x     = 0
    sum_x2    = 0
    sum_v     = 0
    sum_v2    = 0
    sum_rates = 0
    sum_cov   = 0

    !$omp parallel do reduction(+:sum_x,sum_v,sum_v2)
    do ix = 1, num_part
       sum_x = sum_x + pc%particles(ix)%x
       sum_v = sum_v + pc%particles(ix)%v
       sum_v2 = sum_v2 + pc%particles(ix)%v**2
    end do
    !$omp end parallel do

    mean_x = sum_x / num_part
    mean_v = sum_v / num_part
    mean_v2 = sum_v2 / num_part

    !$omp parallel do private(x,v,rates) reduction(+:sum_rates, sum_cov, sum_x2)
    do ix = 1, num_part
       ! Translate to have zero mean
       pc%particles(ix)%x = pc%particles(ix)%x - mean_x

       x = pc%particles(ix)%x
       v = pc%particles(ix)%v
       sum_cov = sum_cov + (x - mean_x) * (v - mean_v)
       sum_x2 = sum_x2 + x**2

       call pc%get_coll_rates(norm2(v), rates)
       sum_rates = sum_rates + rates
    end do
    !$omp end parallel do

    ps%n_samples  = ps%n_samples + num_part
    inv_n_samples = 1.0_dp / ps%n_samples
    fac           = num_part * inv_n_samples

    ! Update bulk velocity
    bulk_v    = mean_x/dt
    ps%bulk_v = ps%bulk_v + (bulk_v - ps%bulk_v) * fac

    ! Update bulk diffusion coefficient
    mean_x2     = sum_x2 / num_part
    bulk_dif    = (mean_x2 - ps%x2_prev) / (2 * dt)
    ps%x2_prev  = mean_x2
    ps%bulk_dif = ps%bulk_dif + (bulk_dif - ps%bulk_dif) * fac

    ! Update flux coefficients
    ps%flux_v  = ps%flux_v + (mean_v - ps%flux_v) * fac
    ps%flux_v2 = ps%flux_v2 + (mean_v2 - ps%flux_v2) * fac
    ! Note that no correction is needed because <x> = 0, see
    ! https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    ps%cov_xv  = ps%cov_xv + sum_cov * corr_fac

    ! Update collision rates
    ps%rates     = ps%rates + (sum_rates - ps%rates*num_part) * inv_n_samples
    ps%coll_rate = sum(ps%rates)
    ps%i_rate    = sum(ps%rates(ionization_colls))
    ps%a_rate    = sum(ps%rates(attachment_colls))

  end subroutine update_particle_stats

  subroutine update_td_from_ps(tds, ps, pc)
    use m_units_constants
    use m_gas
    type(SWARM_td_t), intent(inout) :: tds(:)
    type(part_stats_t), intent(in)  :: ps
    type(PC_t), intent(in)          :: pc
    integer                         :: n
    real(dp)                        :: species_density

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

    do n = 1, size(ix_rates)
       species_density = GAS_get_fraction(trim(pc%cross_secs(n)%gas_name)) * &
            GAS_number_dens
       call update_td(tds(ix_rates(n)), [ps%rates(n) / species_density])
    end do

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
       ! Aim for measurements_per_collision measurements per collision time
       n_measurements = nint(t_measure * ps%coll_rate * &
            measurements_per_collision)
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

       call update_td_from_ps(tds, ps, pc)

       if (verbose > 1) call SWARM_print_results(tds, pc, verbose)

       if (n_swarms >= n_swarms_min) then
          ! Check whether the results are accurate enough
          call get_accuracy(tds, rel_error)
          if (verbose > 0) then
             imax = maxloc(rel_error)
             write(*, '(I6,A40,F12.4)') n_swarms, &
                  trim(tds(imax(1))%description), rel_error(imax(1))
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
    real(dp), intent(out)        :: rel_error(:)
    real(dp)                     :: stddev, mean, var
    integer                      :: i, n
    real(dp), parameter          :: eps = 1e-100_dp

    do i = 1, size(rel_error)
       n = tds(i)%n_measurements
       ! Below we use the norm for vector-based transport data, so that
       ! orientation of the field should not matter
       mean = norm2(tds(i)%val)
       if (n > 1) then
          var = maxval(tds(i)%var / (n-1)) / n ! Maximal variance
       else
          var = maxval(tds(i)%var / (n)) / n ! Maximal variance
       end if
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

  subroutine SWARM_print_results(tds, pc, verbose)
    use m_units_constants
    use m_gas
    type(SWARM_td_t), intent(in) :: tds(:) !< The transport data
    type(pc_t), intent(in)       :: pc
    integer, intent(in)          :: verbose
    integer                      :: i, i_dim
    real(dp)                     :: fac, tmp, std, N0, val
    real(dp)                     :: energy, mu, rel_error(size(tds))
    logical                      :: magnetic_field_used
    character(len=2)             :: dimnames(3) = [" x", " y", " z"]

    N0 = GAS_number_dens
    i = tds(1)%n_measurements

    if (i > 1) then
       fac = sqrt(1.0_dp / (i * (i-1)))
    else
       fac = 1.0_dp             ! To prevent division by zero
    end if

    magnetic_field_used = (abs(SWARM_field%Bz) > 0.0_dp)

    call print_td("Gas number density (1/m3)", N0, 0.0_dp)
    call print_td("Electric field (V/m)", &
         norm2(SWARM_field%E_vec), 0.0_dp)
    call print_td("Electric field / N (Td)", &
         norm2(SWARM_field%E_vec)/N0*1e21_dp, 0.0_dp)
    if (magnetic_field_used) then
       call print_td("Magnetic field (T)", SWARM_field%Bz, 0.0_dp)
       call print_td("E-B field angle (degree)", &
            SWARM_field%angle_deg, 0.0_dp)
       call print_td("omega_c (radian/s)", SWARM_field%omega_c, 0.0_dp)
    end if

    ! mean energy
    tmp    = 0.5_dp * UC_elec_mass / UC_elec_volt
    energy = tmp * sum(tds(ix_flux_v2)%val)
    std    = tmp * sqrt(sum(tds(ix_flux_v2)%var)) * fac
    call print_td("Mean energy (eV)", energy, std)

    ! Ionization coefficient
    val = tds(ix_alpha)%val(1)/N0
    std = fac * sqrt(tds(ix_alpha)%var(1))/N0
    call print_td("Townsend ioniz. coef. alpha/N (m2)", val, std)

    ! Attachment coefficient
    val = tds(ix_eta)%val(1)/N0
    std = fac * sqrt(tds(ix_eta)%var(1))/N0
    call print_td("Townsend attach. coef. eta/N (m2)", val, std)

    if (magnetic_field_used) then
       ! mobility parallel to E
       mu = -dot_product(tds(ix_flux_v)%val,  SWARM_field%E_vec) &
            / max(epsilon(1.0_dp), sum(SWARM_field%E_vec**2))
       std = fac * sqrt(dot_product(tds(ix_flux_v)%var,  SWARM_field%E_vec)) &
            / max(epsilon(1.0_dp), norm2(SWARM_field%E_vec)**1.5_dp)
       call print_td("mu_E (m2/V/s)", mu, std)

       ! mobility parallel to B
       if (abs(SWARM_field%Ez) > sqrt(epsilon(1.0_dp))) then
          mu = -tds(ix_flux_v)%val(3) / SWARM_field%Ez
          std = fac * sqrt(tds(ix_flux_v)%var(3)) / abs(SWARM_field%Ez)
       else
          mu = 0
          std = 0
       end if
       call print_td("mu_B (m2/V/s)", mu, std)

       ! mobility perpendicular to B (y-velocity over Ey)
       if (abs(SWARM_field%Ey) > sqrt(epsilon(1.0_dp))) then
          mu = -tds(ix_flux_v)%val(2) / SWARM_field%Ey
          std = fac * sqrt(tds(ix_flux_v)%var(2)) / abs(SWARM_field%Ey)
       else
          mu = 0
          std = 0
       end if
       call print_td("mu_xB (m2/V/s)", mu, std)

       ! mobility in ExB-direction (x-velocity over Ey = E_perp)
       if (abs(SWARM_field%Ey) > sqrt(epsilon(1.0_dp)) .and. &
            abs(SWARM_field%Bz) > sqrt(epsilon(1.0_dp))) then
          mu = tds(ix_flux_v)%val(1) / SWARM_field%Ey
          std = fac * sqrt(tds(ix_flux_v)%var(1)) / abs(SWARM_field%Ey)
       else
          mu = 0
          std = 0
       end if
       call print_td("mu_ExB (m2/V/s)", mu, std)
    else
       ! mobility parallel to E
       mu = -dot_product(tds(ix_flux_v)%val,  SWARM_field%E_vec) &
            / max(epsilon(1.0_dp), sum(SWARM_field%E_vec**2))
       std = fac * sqrt(dot_product(tds(ix_flux_v)%var,  SWARM_field%E_vec)) &
            / max(epsilon(1.0_dp), norm2(SWARM_field%E_vec)**1.5_dp)
       call print_td("Mobility (m2/V/s)", mu, std)
       call print_td("Mobility *N (1/m/V/s)", mu*N0, std*N0)

       mu = -dot_product(tds(ix_bulk_v)%val,  SWARM_field%E_vec) &
            / max(epsilon(1.0_dp), sum(SWARM_field%E_vec**2))
       std = fac * sqrt(dot_product(tds(ix_bulk_v)%var,  SWARM_field%E_vec)) &
            / max(epsilon(1.0_dp), norm2(SWARM_field%E_vec)**1.5_dp)
       call print_td("Bulk mobility (m2/V/s)", mu, std)
       call print_td("Bulk mobility *N (1/m/V/s)", mu*N0, std*N0)

       ! Longitudinal diffusion
       val = tds(ix_flux_diff)%val(3)
       std = sqrt(tds(ix_flux_diff)%var(3))
       call print_td("Flux L diffusion coef. (m2/s)", val, std)
       call print_td("Flux L diffusion coef. *N (1/m/s)", val*N0, std*N0)

       val = tds(ix_bulk_diff)%val(3)
       std = sqrt(tds(ix_bulk_diff)%var(3))
       call print_td("Bulk L diffusion coef. (m2/s)", val, std)
       call print_td("Bulk L diffusion coef. *N (1/m/s)", val*N0, std*N0)

       ! Transverse diffusion
       val = 0.5_dp * sum(tds(ix_flux_diff)%val(1:2))
       std = sqrt(0.5_dp * sum(tds(ix_flux_diff)%var(1:2)))
       call print_td("Flux T diffusion coef. (m2/s)", val, std)
       call print_td("Flux T diffusion coef. *N (1/m/s)", val*N0, std*N0)

       val = 0.5_dp * sum(tds(ix_bulk_diff)%val(1:2))
       std = sqrt(0.5_dp * sum(tds(ix_bulk_diff)%var(1:2)))
       call print_td("Bulk T diffusion coef. (m2/s)", val, std)
       call print_td("Bulk T diffusion coef. *N (1/m/s)", val*N0, std*N0)
    end if

    call get_accuracy(tds, rel_error)

    ! The other swarm parameters
    do i = 1, size(tds)
       if (tds(i)%n_dim == 1) then
          call print_td(trim(tds(i)%description) // " (" // &
               trim(tds(i)%unit) // ")", tds(i)%val(1), &
               sqrt(tds(i)%var(1)) * fac, rel_error(i))
       else
          do i_dim = 1, tds(i)%n_dim
             call print_td(trim(tds(i)%description) // dimnames(i_dim) &
                  // " (" // trim(tds(i)%unit) // ")", &
                  tds(i)%val(i_dim), sqrt(tds(i)%var(i_dim)) * fac, &
                  rel_error(i))
          end do
       end if
    end do

  end subroutine SWARM_print_results

  subroutine print_td(name, val, stddev, convergence)
    character(len=*), intent(in)   :: name
    real(dp), intent(in)           :: val, stddev
    real(dp), intent(in), optional :: convergence
    character(len=40)              :: name_left
    character(len=*), parameter    :: fmtstr = "(A40,ES14.6,ES10.2,F7.2)"

    name_left = name
    if (present(convergence)) then
       write(*, fmtstr) name_left, val, stddev, convergence
    else
       write(*, fmtstr) name_left, val, stddev, 0.0_dp
    end if
  end subroutine print_td

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
