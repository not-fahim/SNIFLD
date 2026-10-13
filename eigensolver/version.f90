module version_mod
        implicit none

contains
        subroutine print_version_data(date_str, time_str)
                character(len=8), intent(out)  :: date_str
                character(len=8), intent(out)  :: time_str
        
                print *, "Code name             : SNIF Sn in Fortran"
                print *, "Author                : Md. Fahim Foraejee"
                print *, "Version               : 0.0.1"
                call date_and_time(date = date_str, &
                                   time = time_str)
                print *, "Execution  date       : ", date_str(1:4), ":", date_str(5:6), ":", date_str(7:8)
                print *, "Time                  : ", time_str(1:2), ":", time_str(3:4), ":", time_str(5:6) 

         end subroutine print_version_data
         
  subroutine write_version_to_file(input_name, date_str, time_str)
    character(len=8), intent(in)  :: date_str
    character(len=8), intent(in)  :: time_str
    character(len = 255), intent(in) :: input_name 
    !  ---local
    character(len = 259) :: output_name
    integer :: out_id, out_io_stat


     output_name = trim(input_name) // ".out"
     out_id = 71


    ! --- open output it replaces if already eixts ---
    open(unit = out_id, file = output_name, action = 'write',  status = 'replace', iostat = out_io_stat)


    if(out_io_stat /= 0) then
      print *, "Error: Could not open output file: ", output_name
    stop
    end if
    write(out_id, *) "Code name             : SNIF Sn in Fortran"
    write(out_id, *) "Author                : Md. Fahim Foraejee"
    write(out_id, *) "Version               : 0.0.1"
    write(out_id, *) "Execution  date       : ", date_str(1:4), ":", date_str(5:6), ":", date_str(7:8)
    write(out_id, *) "Time                  : ", time_str(1:2), ":", time_str(3:4), ":", time_str(5:6) 
    close(out_id)
  end subroutine write_version_to_file
        
end module version_mod 
