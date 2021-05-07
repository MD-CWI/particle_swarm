program test_thermalization
  use m_particle_core
  use m_cross_sec
  use m_units_constants

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  character(len=*), parameter :: cs_file = "cs_constant_elastic.txt"
  character(len=*), parameter :: gas_name = "N2"

  integer, parameter      :: max_num_part  = 1000*1000
  integer, parameter      :: init_num_part = 1000
  integer, parameter      :: max_num_steps = 100
  integer, parameter      :: lkp_tbl_size  = 1000
  real(dp), parameter     :: delta_t       = 5.0e-9_dp
  real(dp), parameter     :: max_en_eV     = 10.0_dp
  real(dp), parameter     :: neutral_dens  = 2.5e25_dp
  real(dp), parameter     :: part_mass     = UC_elec_mass
  real(dp), parameter     :: init_accel(3) = (/0.0_dp, 0.0_dp, 0.0_dp/)
  real(dp)                :: pos(3), vel(3), accel(3), weight
  integer                 :: ll, step, num_colls
  type(CS_t), allocatable :: cross_secs(:)
  type(PC_t)              :: pc

  ! Use constant momentum transfer cross section
  print *, "Reading in cross sections from ", trim(cs_file)
  call CS_add_from_file(cs_file, gas_name, neutral_dens, max_en_eV, &
       cross_secs)
  call pc%initialize(part_mass, max_num_part)
  call pc%use_cross_secs(max_en_eV, lkp_tbl_size, cross_secs)

  pc%gas_temperature = 300.0_dp

  num_colls = pc%get_num_colls()
  deallocate(cross_secs)

  print *, "Creating initial particles"
  do ll = 1, init_num_part
     pos    = 0.0_dp
     call random_number(vel)
     vel    = (vel - 0.5_dp) * 1e5
     accel  = init_accel
     weight = 1
     call pc%create_part(pos, vel, accel, weight, 0.0_dp)
  end do

  print *, "step time        mean_eV     mean/expected"
  do step = 1, max_num_steps
     call print_stats(step, (step-1) * delta_t)
     call pc%advance_openmp(delta_t)
  end do

  call print_stats(step, (step-1) * delta_t)

contains

  subroutine part_stats(part, vec)
    type(PC_part_t), intent(in) :: part
    real(dp), intent(out) :: vec(:)
    vec(1:3) = part%w * part%x
    vec(4:6) = part%w * part%v
    vec(7:9) = part%w * part%a
    vec(10) = part%w * PC_v_to_en(part%v, part_mass)
    vec(11) = part%w
  end subroutine part_stats

  subroutine print_stats(step, time)
    integer, intent(in)  :: step
    real(dp), intent(in) :: time
    integer              :: n_part
    real(dp)             :: sum_x(3), sum_v(3), sum_a(3), sum_en
    real(dp)             :: sum_weight, mean_eV, expected_eV
    real(dp)             :: sum_vec(11)

    n_part = 0

    n_part = n_part + pc%get_num_sim_part()
    call pc%compute_vector_sum(part_stats, sum_vec)

    sum_weight = sum_vec(11)
    sum_x = sum_vec(1:3)
    sum_v = sum_vec(4:6)
    sum_a = sum_vec(7:9)
    sum_en =  sum_vec(10)

    mean_eV = sum_en / (sum_weight * UC_elec_volt)
    expected_eV = UC_boltzmann_const * pc%gas_temperature * 1.5_dp / UC_elec_volt
    write(*, "(I4, E12.4, E12.4, F8.4)") step, time, &
         mean_eV, mean_eV/expected_eV
  end subroutine print_stats

end program test_thermalization
