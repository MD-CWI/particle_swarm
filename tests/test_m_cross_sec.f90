program test_m_cross_sec
   use m_cross_sec

   integer, parameter :: dp = kind(0.0d0)
   type(CS_t), allocatable :: cross_secs(:)
   ! Test 1 variables
   type(CS_t), allocatable :: TEST_1_cross_secs(:)
   real(dp)                    :: TEST_1_req_energy = 1.0e5_dp
   ! Test 2 variables
   type(CS_t), allocatable :: TEST_2_cross_secs(:)
   real(dp)                    :: TEST_2_req_energy = 10100
   ! Test 3 variables
   type(CS_t), allocatable :: TEST_3_cross_secs(:)
   real(dp)                    :: TEST_3_req_energy = 1.0e5_dp
   ! Test 4 variables
   type(CS_t), allocatable :: TEST_4_cross_secs(:)
   ! Test 5 variables
   type(CS_t), allocatable :: TEST_5_cross_secs(:)
   ! Test 6 variables
   type(CS_t), allocatable :: TEST_6_cross_secs(:)



   character(len=*), parameter :: first_in = "test_m_cross_sec_input.txt"
   character(len=*), parameter :: first_out = "test_m_cross_sec_all_output.txt"
   character(len=*), parameter :: first_summ = "test_m_cross_sec_summary_output.txt"
   character(len=*), parameter :: second_out = "test_m_cross_sec_all_2_output.txt"
   character(len=*), parameter :: gas_1 = "Ar", gas_2 = "N2", gas_3 = "O2"

   call CS_add_from_file(first_in, gas_1, 1.0_dp, 1.0e3_dp, cross_secs, &
        handle_effective=CS_effective_keep)
   call CS_add_from_file(first_in, gas_2, 1.0_dp, 1.0e3_dp, cross_secs, &
        handle_effective=CS_effective_keep)
   call CS_add_from_file(first_in, gas_3, 1.0_dp, 1.0e3_dp, cross_secs, &
        handle_effective=CS_effective_keep)

   call CS_write_summary(cross_secs, first_summ)
   print *, "First pass: read in ", size(cross_secs), " cross sections"

   print *, "Out of bounds input data tests."
   ! Throw error when called, but there does not seem to be nice exception handling so cannot run this in test
   !call CS_add_from_file(first_in, gas_1, 1.0_dp, 1.0e5_dp, cross_secs, opt_out_of_bounds_upper=CS_extrapolate_Throw)

   print *, "TEST 1: Out of bounds upper: CS_extrapolate_zero: If input data energy does not reach required energy add a 0 cs line."
   call CS_add_from_file(first_in, gas_1, 1.0_dp, TEST_1_req_energy, &
        TEST_1_cross_secs, opt_out_of_bounds_upper=CS_extrapolate_zero, &
        handle_effective=CS_effective_keep)
   print *, "Required energy: ", TEST_1_req_energy
   print *, "Last energy: ", TEST_1_cross_secs(1)%en_cs(1, TEST_1_cross_secs(1)%n_rows), " expected to be: ",&
                             TEST_1_cross_secs(1)%en_cs(1, TEST_1_cross_secs(1)%n_rows -1 ) * (1 + epsilon(1.0_dp))
   print *, "Last cross section: ", TEST_1_cross_secs(1)%en_cs(2, TEST_1_cross_secs(1)%n_rows), " expected to be: 0"

   print *, "TEST 2: Out of bounds upper: CS_extrapolate_linear: If input data energy does not reach required energy ", &
            "linearly interpolate to it."
   call CS_add_from_file(first_in, gas_1, 1.0_dp, TEST_2_req_energy, &
        TEST_2_cross_secs, opt_out_of_bounds_upper=CS_extrapolate_linear, &
        handle_effective=CS_effective_keep)
   print *, "Required energy: ", TEST_2_req_energy
   print *, "Last energy: ", TEST_2_cross_secs(1)%en_cs(1, TEST_2_cross_secs(1)%n_rows), " expected to be: 10100"
   print *, "Last cross section: ", TEST_2_cross_secs(1)%en_cs(2, TEST_2_cross_secs(1)%n_rows), &
             " expected to be: 1.7316666666666668e-21"

   print *, "TEST 3: Out of bounds upper: CS_extrapolate_linear: If input data energy does not reach required energy ", &
            "linearly interpolate to it BUT if interpolated value is negative warn user and set cs to 0."
   call CS_add_from_file(first_in, gas_1, 1.0_dp, TEST_3_req_energy, &
        TEST_3_cross_secs, opt_out_of_bounds_upper=CS_extrapolate_linear, &
        handle_effective=CS_effective_keep)
   print *, "Required energy: ", TEST_3_req_energy
   print *, "Last energy: ", TEST_3_cross_secs(1)%en_cs(1, TEST_3_cross_secs(1)%n_rows), " expected to be: 100000"
   print *, "Last cross section: ", TEST_3_cross_secs(1)%en_cs(2, TEST_3_cross_secs(1)%n_rows), " expected to be: 0"

   print *, "TEST 4: Out of bounds lower: CS_extrapolate_constant: Keep the first (energy, cs) pair as it is."
   call CS_add_from_file(first_in, gas_1, 1.0_dp, 100.0_dp, TEST_4_cross_secs, &
        opt_out_of_bounds_lower=CS_extrapolate_constant, &
        handle_effective=CS_effective_keep)
   print *, "First energy: ", TEST_4_cross_secs(1)%en_cs(1, 1), " expected to be: 0"
   print *, "First cross section: ", TEST_4_cross_secs(1)%en_cs(2, 1), " expected to be: 7.5000e-20"

   print *, "TEST 5: Out of bounds lower: CS_extrapolate_zero: Add an extra line at the start of the input data ", &
            "with energy 1% less than first point in the input data and the cross section 0"
   call CS_add_from_file(first_in, gas_1, 1.0_dp, 100.0_dp, TEST_5_cross_secs, &
        opt_out_of_bounds_lower=CS_extrapolate_zero, &
        handle_effective=CS_effective_keep)
   print *, "First energy: ", TEST_5_cross_secs(4)%en_cs(1, 1), " expected to be: ", &
                              TEST_5_cross_secs(4)%en_cs(1, 2) * (1 - epsilon(1.0_dp))
   print *, "First cross section: ", TEST_5_cross_secs(4)%en_cs(2, 1), " expected to be: 0"

   print *, "TEST 6: Out of bounds lower: CS_extrapolate_zero: SHOULD NOT Add an extra line at the start of the input data ", &
            "with energy 1% less than first point in the input data and the cross section 0"
   call CS_add_from_file(first_in, gas_1, 1.0_dp, 100.0_dp, TEST_6_cross_secs, &
        opt_out_of_bounds_lower=CS_extrapolate_zero, &
        handle_effective=CS_effective_keep)
   print *, "First energy: ", TEST_6_cross_secs(1)%en_cs(1, 1), " expected to be: 0"
   print *, "First cross section: ", TEST_6_cross_secs(1)%en_cs(2, 1), " expected to be: 7.5000e-20"


   deallocate(cross_secs)
end program test_m_cross_sec
