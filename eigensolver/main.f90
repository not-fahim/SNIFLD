program solver_app
      use global_vars
      use parser_mod           
      use check_mod
      use version_mod
      use output_mod
      use transport_mod
!     use ddsolver_mod
!     use sweep_mod
      implicit none

      integer :: n_args
      call cpu_time(start)

      
      n_args = command_argument_count()

      if (n_args /= 1) then
         print *, "Usage: ./solver_app <file>"
         stop
      end if
      
      call get_command_argument(1, input_name)
      
      call print_version_data(date_str, time_str)
      call write_version_to_file(input_name, date_str, time_str)

      ! the parser will read the file, allocate these arrays fill them with data.
      call read_input_file(input_name)
      ! check physicality of the input

      call check_input()
      ! write input data to file

      call  write_input_to_file (input_name) 
!      print *, "will solve D.O. here"
      call transport_solver() 
      call write_phi_to_file(input_name, start)
     
      
end program solver_app
