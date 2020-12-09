module m_io

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)

   public :: IO_write_td
   public :: IO_write_td_cols
   public :: IO_write_coeffs

contains

   subroutine IO_write_td(tData, tNames, gasName, filename, EfieldIx, energyIx)
      real(dp), intent(IN)         :: tData(:, :)
      character(LEN=*), intent(IN) :: tNames(:), gasName, filename
      integer, intent(IN)          :: EfieldIx, energyIx

      integer :: n, nE, i_max_td, i_max_fld
      integer :: ios, OutUnit = 456

      open(UNIT=OutUnit, FILE=filename, IOSTAT=ios)

      i_max_fld = ubound(tData, 2)
      i_max_td  = ubound(tData, 1)

      ! Write data depending on Efield
      do n = 1, i_max_td
         if (n == EfieldIx) cycle
         write(OutUnit, *) ""
         write(OutUnit, *) trim(tNames(EfieldIx))//"_vs_"//trim(tNames(n))
         write(OutUnit, *) trim(gasName)
         write(OutUnit, *) "------------------------"
         do nE = 1, i_max_fld
            write(OutUnit, *) tData(EfieldIx, nE), tData(n, nE)
         end do
         write(OutUnit, *) "------------------------"
      end do

      ! Write data depending on energy
      do n = 1, i_max_td
         if (n == energyIx) cycle
         write(OutUnit, *) ""
         write(OutUnit, *) trim(tNames(energyIx))//"_vs_"//trim(tNames(n))
         write(OutUnit, *) trim(gasName)
         write(OutUnit, *) "------------------------"
         do nE = 1, i_max_fld
            write(OutUnit, *) tData(energyIx, nE), tData(n, nE)
         end do
         write(OutUnit, *) "------------------------"
      end do

      write(OutUnit, *) ""
      close(UNIT=OutUnit)
   end subroutine IO_write_td

   ! Write transport date into a file which can be plotted more easily
   subroutine IO_write_td_cols(tData, tNames, gasName, filename)
      real(dp), optional,intent(IN) :: tData(:, :)
      character(LEN=*), optional,intent(IN) :: tNames(:), gasName, filename

      integer :: n, nE, i_max_td, i_max_fld

      integer :: ios, OutUnit = 457

      open(UNIT=OutUnit, FILE=filename, IOSTAT=ios)

      i_max_fld = ubound(tData, 2)
      i_max_td  = ubound(tData, 1)

      ! Write data depending on Efield
      write(OutUnit, *) "# Swarm data"
      write(OutUnit, *) "# Gas:" // trim(gasName)
      do n=1, i_max_td
         write(OutUnit, ADVANCE="NO", FMT="(A)") " #" // trim(tNames(n))
      end do
      write(OutUnit, *) "#"
      write(OutUnit, *) "# ------------------------"
      do nE = 1, i_max_fld
         write(OutUnit, *) tData(:, nE)
      end do

   end subroutine IO_write_td_cols

   subroutine IO_write_coeffs(pc, filename)
      use m_particle_core
      type(PC_t), intent(in) :: pc
      integer                              :: n_coeffs
      character(len=*), intent(in)         :: filename

      real(dp), allocatable                :: coeff_data(:,:)
      character(len=20), allocatable :: coeff_names(:)

      call pc%get_coeffs(coeff_data, coeff_names, n_coeffs)

      if (n_coeffs > 0) then
         call write_data_2d(filename, coeff_data, coeff_names, 30, do_transpose = .true.)
      end if
   end subroutine IO_write_coeffs

   subroutine write_data_2d(filename, data_2d, col_names, col_width, do_transpose)
      character(len=*), intent(in)        :: filename
      real(dp), intent(in)                :: data_2d(:,:)
      integer, intent(in)                 :: col_width
      character(len=*), intent(in) :: col_names(:)
      logical, intent(in), optional       :: do_transpose

      integer                             :: my_unit, n, n_rows, n_cols, io_state
      character(len=20)             :: fmt_string
      real(dp), allocatable               :: copy_of_data(:,:)
      logical                             :: transpose_data

      my_unit = 333
      if (present(do_transpose)) then
         transpose_data = do_transpose
      else
         transpose_data = .false.
      end if

      if (transpose_data) then
         n_rows = size(data_2d, 2)
         n_cols = size(data_2d, 1)
         allocate(copy_of_data(n_rows, n_cols))
         copy_of_data = transpose(data_2d)
      else
         n_rows = size(data_2d, 1)
         n_cols = size(data_2d, 2)
         allocate(copy_of_data(n_rows, n_cols))
         copy_of_data = data_2d
      end if

      if (size(col_names) /= n_cols) then
         print *, "write_data_2d: incompatible argument sizes"
         stop
      end if

      open(unit=my_unit, file=filename, iostat=io_state)
      if (io_state /= 0) then
         print *, "write_data_2d: error writing " // filename
         stop
      end if

      ! Write header
      do n = 1, n_cols
         write(my_unit, advance="NO", fmt="(A)") '"' // &
              trim(col_names(n)) // '" '
      end do

      write(my_unit, *) ""

      ! Create format string for data
      write(fmt_string, fmt="(A,I0,A,I0,A,I0,A)") "(", n_cols, "E", col_width, ".", col_width - 9, "e3)"

      ! Write data
      do n = 1, n_rows
         write(my_unit, fmt_string) copy_of_data(n, :)
      end do

      close(my_unit)

   end subroutine write_data_2d

end module m_io
