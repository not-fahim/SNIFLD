module parser_mod
      use global_vars
      implicit none
contains

      subroutine read_input_file(filename)

        character(len = 255), intent(in)        :: filename
        integer                                 :: io_stat       ! status of IO operations
        integer                                 :: file_id
        

        ! --- internal parsing variables ---
        character(len=1000) :: line_buf   ! buffer to read text lines
        integer             :: vals_on_line ! how many numbers on the line
        integer             :: n    ! loop counters
        integer             :: current_idx  ! tracking index for vectors
        integer             :: read_count! tracking number of elements read

        !  open File
        open(newunit=file_id, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "Error: Could not open file ", trim(filename)
            stop
        end if

        print *, "Reading file: ", trim(filename)

        !  read grid dimensions (Nx, Ny) .........................
        vals_on_line = 0
        do while(vals_on_line == 0)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
        end do

        if(vals_on_line /= 1) then
                print *, "Error in card 1: Grid dimension not provided correctly"
                print *, "expected: 1 values: Nx , but found", vals_on_line, "values"
                stop
        end if

        read(line_buf, *) Nx
       ! print *, "Grid:", I, "x", J

        !  Read delx ...............................................
        allocate(delx(Nx))
        current_idx = 1
        
        do while (current_idx <= Nx)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
           
           ! empty line  
           if (vals_on_line == 0) cycle

           ! did not follow input rule
           if(vals_on_line > 5) then
                print *, "Error in card 2: in delx no more than 5 vals are allowed per line"
                stop
           end if
           ! upto what idx to write 
           read_count = current_idx + vals_on_line - 1
           if(read_count > Nx) then
                print *, "Error in card 2: Too many vals for cell x dimension"
                print *, "Expected:", Nx , " found:", read_count          
                stop
           end if
 

           read(line_buf, *) delx(current_idx : read_count)
           current_idx = read_count + 1 
        end do

        if(read_count /= Nx) then

              print *, "Error in card 2: Expected", Nx, "cell x dimensions, found:", read_count
              stop
        end if

        ! read quadrature data..............................................
        vals_on_line = 0
        do while(vals_on_line == 0)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
        end do
        if(vals_on_line /= 1) then
                print *, "Error in card 3: nnumber of angles/octant not provided correctly"
                print *,  "expected:1 value: Nq, but found", vals_on_line, "values" 
                stop
        end if


        read(line_buf, *) Nq

        allocate(mu(2*Nq))
        allocate(w(2*Nq))

        n = 1 
        do while(n<=Nq)  
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
           
           ! empty line  
           if (vals_on_line == 0) cycle

           ! did not follow input rule
           if(vals_on_line /= 2) then
                print *, "Error in card 4: in ordinates each direction should have 2 reals per line"
                stop
           end if
           
           read(line_buf, *) mu(n), w(n) 
           n = n + 1
        end do

        do n = 1, Nq 
            mu(Nq+n) = -mu(n)
            w(Nq+n) = w(n)
        end do

        Nq = 2*Nq


        ! read materials (Nm, then sigmas) ....................................
        vals_on_line = 0
        do while(vals_on_line == 0)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
        end do
        if(vals_on_line /= 1) then
                print *, "Error in card 5: number of material regions are not provided correctly"
                print *,  "expected:1 value Nm, but found", vals_on_line, "values" 
                stop
        end if


        read(line_buf, *) Nm

        allocate(sigma_t(Nm), sigma_s(Nm), nu(Nm), sigma_f(Nm))
        n = 1
        do while( n <= Nm)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
           
           ! empty line  
           if (vals_on_line == 0) cycle

           ! did not follow input rule
           if(vals_on_line /= 4) then
                print *, "Error in card 6: each material should have 4 reals per line"
                print *, "expected format: sigmaT sigmaS nu sigmaF; found:", vals_on_line
                stop
           end if
 
            read(line_buf, *) sigma_t(n), sigma_s(n), nu(n), sigma_f(n)
            n = n+1
        end do

        ! BC left ............................................. 
        Nbc = Nq/2
        allocate(left_bc(Nbc))
        current_idx = 1
        do while (current_idx <= Nbc)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
        
           ! empty line  
           if (vals_on_line == 0) cycle

           ! did not follow input rule
           if(vals_on_line > 10) then
               print *, "Error in card 7.1, left BC: no more than 10 reals are  allowed per line"
               stop
           end if

           ! upto what idx to write 
           read_count = current_idx + vals_on_line - 1
           if(read_count > Nbc) then
               print *, "Error in card 7.1, left BC: Too many args for left BC data"
               print *, "Expected:", Nbc , " found:", read_count          
               stop
           end if
 
           read(line_buf, *) left_bc(current_idx : read_count)
           current_idx = read_count + 1 
        end do

        ! BC right ............................................. 
        allocate(right_bc(Nbc))
        current_idx = 1
        do while (current_idx <= Nbc)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
        
           ! empty line  
           if (vals_on_line == 0) cycle

           ! did not follow input rule
           if(vals_on_line > 10) then
               print *, "Error in card 7.2, right BC: no more than 10 reals are  allowed per line"
               stop
           end if

           ! upto what idx to write 
           read_count = current_idx + vals_on_line - 1
           if(read_count > Nbc) then
               print *, "Error in card 7.2, right BC: Too many args for right BC data"
               print *, "Expected:", Nbc , " found:", read_count          
               stop
           end if
 
           read(line_buf, *) right_bc(current_idx : read_count)
           current_idx = read_count + 1 
        end do


        ! material map .................................
        allocate(matmap(Nx))

        current_idx = 1
        
        do while (current_idx <= Nx)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
        
           ! empty line  
           if (vals_on_line == 0) cycle

           ! did not follow input rule
           if(vals_on_line > 10) then
               print *, "Error in card 8: in matmap no more than 10 reals are  allowed per line"
               stop
           end if

           ! upto what idx to write 
           read_count = current_idx + vals_on_line - 1
           if(read_count > Nx) then
               print *, "Error in card 8: Too many args for material map data"
               print *, "Expected:", Nx , " found:", read_count          
               stop
           end if
 
           read(line_buf, *) matmap(current_idx : read_count)
           current_idx = read_count + 1 
        end do



        ! source zeorth moment map ..............................
        allocate(q_ext_i(Nx))
        current_idx = 1
        
        do while (current_idx <= Nx)
          call read_clean_line(file_id, line_buf)
          vals_on_line = count_items(line_buf)
        
          ! empty line  
          if (vals_on_line == 0) cycle

          ! did not follow input rule
          if(vals_on_line > 10) then
              print *, "Error in card 9: in source map no more than 10 vals are allowed per line"
              stop
          end if
          ! upto what idx to write 
          read_count = current_idx + vals_on_line - 1
          if(read_count > Nx) then
              print *, "Error in card 9: Too many args for src map"
              print *, "Expected:", Nx , " found:", read_count          
              stop
          end if
 
          read(line_buf, *) q_ext_i(current_idx : read_count)
          current_idx = read_count + 1 
        end do


        ! source first spatial  moment map ..............................
        allocate(q_ext_hat_i(Nx))
        current_idx = 1
        
        do while (current_idx <= Nx)
          call read_clean_line(file_id, line_buf)
          vals_on_line = count_items(line_buf)
        
          ! empty line  
          if (vals_on_line == 0) cycle

          ! did not follow input rule
          if(vals_on_line > 10) then
              print *, "Error in card 9.1: in source map no more than 10 vals are allowed per line"
              stop
          end if
          ! upto what idx to write 
          read_count = current_idx + vals_on_line - 1
          if(read_count > Nx) then
              print *, "Error in card 9.1: Too many args for src map"
              print *, "Expected:", Nx , " found:", read_count          
              stop
          end if
 
          read(line_buf, *) q_ext_hat_i(current_idx : read_count)
          current_idx = read_count + 1 
        end do

        !  read inner tolerance and max inner iterations (inner_tol, max_inners) .........................
        vals_on_line = 0
        do while(vals_on_line == 0)
           call read_clean_line(file_id, line_buf)
           vals_on_line = count_items(line_buf)
        end do

        if(vals_on_line /= 2) then
                print *, "Error in card 10: in parsing convergence and stopping criteria"
                print *, "expected: 2 values: convergence criterion, stopping criterion  but found", vals_on_line, "values"
                stop
        end if

        read(line_buf, *) inner_tol, max_inners


        close(file_id)
        print *, "Input file read successfully."


      end subroutine read_input_file




      ! HELPER: reads a line and strips '#' comments 
      subroutine read_clean_line(fid, buffer)
          integer, intent(in) :: fid
          character(len=*), intent(out) :: buffer
          integer :: io, comment_pos

          read(fid, '(A)', iostat=io) buffer
          if (io /= 0) then
              print *, "Error: Unexpected end of file."
              print *, "This error results when fewer arguments are provided for an array"
              print *, "Check matmap, srcmap, src moment map"
              stop
          end if

          ! strip comments
          comment_pos = index(buffer, '#')
          if (comment_pos > 0) buffer(comment_pos:) = ' '
      end subroutine read_clean_line




      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! HELPER: counts number of values in line 
      integer function count_items(s)
         character(len=*), intent(in) :: s
         integer                      :: k
         logical                      :: in_number
         
         count_items = 0
         in_number = .false.
         
         do k = 1, len_trim(s)
            if (s(k:k) /= ' ' .and. s(k:k) /= char(9)) then !char(9) is tab
                if (.not. in_number) then
                    count_items = count_items + 1
                    in_number = .true.
                end if
            else
                in_number = .false.
            end if
         end do
      end function count_items

end module parser_mod
