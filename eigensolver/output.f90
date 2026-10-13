module output_mod
    use global_vars
    implicit none

contains
    subroutine write_input_to_file(input_name)

        
        character(len = 255), intent(in) :: input_name 
        !  ---local
        character(len = 259) :: output_name
        integer :: out_id, out_io_stat
        integer ::  n, ii, mat_id


        output_name = trim(input_name) // ".out"
        out_id = 51


        ! --- open output it replaces if already eixts ---

        open(unit = out_id, file = output_name, action = 'write', &
             status = 'unknown', position = 'append', iostat = out_io_stat)


        if(out_io_stat /= 0) then
             print *, "Error while printing inputs: Could not open output file: ", output_name
             stop

        end if
        
        ! angles per octant
        write(out_id, *)

        write(out_id, '(A)') "Discrete Ordinates/Octant"

        write(out_id, '(A)') "    n            mu             w"
        
        do n = 1, Nq
            write(out_id, '(I5, 2(2x, ES16.8))') n, mu(n), w(n)
        end do

        ! boundary
        write(out_id, *) ! Empty line for spacing

        write(out_id, '(A)') "Boundary Conditions"
        write(out_id, '(A)') "Left"
        write(out_id, '(A)') "    n            mu             psi_m"
        
        do n=1,Nbc
            write(out_id, '(I5, 2(2x, ES16.8))') n, mu(n), left_bc(n)
        end do

        write(out_id, '(A)') "Right"
        write(out_id, '(A)') "    n            mu             psi_m"
        
        do n=1,Nbc
            write(out_id, '(I5, 2(2x, ES16.8))') n, mu(n+Nbc), right_bc(n)
        end do
          
        

        ! computational cell data
        write(out_id, *) !skip line
        write(out_id, '(A)') "Computational Cell Data"
        ! Printing a descriptive header to match column widths
        write(out_id, '(A)') "    i  Mat    dx         SigmaT          SigmaS         nu        SigmaF          q          q_hat"
  

        do ii = 1, Nx
            mat_id = matmap(ii)
            
            ! Notes on Format:
            ! 2I6: Three integers (i,  mat), width 6 each
            ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
            write(out_id, '(I5,I3, 7(1x, ES13.6))') &
                ii, mat_id, &
                delx(ii), &
                sigma_t(mat_id), sigma_s(mat_id), &
                nu(mat_id), sigma_f(mat_id), &
                q_ext_i(ii), &
                q_ext_hat_i(ii)
        end do
        write(out_id, *) !skip line
        write(out_id, '(A)') "will solve D.O. here"
        close(out_id)

    end subroutine write_input_to_file
    
subroutine write_phi_to_file(input_name, start)

   character(len = 255), intent(in) :: input_name 
   !  ---local
   character(len = 259) :: output_name
   integer :: out_id, out_io_stat
   integer ::  ii, mm
   real    :: exec_time, start
   real(kind(0d0)) :: x_val

   output_name = trim(input_name) // ".out"
   out_id = 61


   ! --- open output it replaces if already eixts ---

   open(unit = out_id, file = output_name, action = 'write', &
     status = 'old', position = 'append',  iostat = out_io_stat)


   if(out_io_stat /= 0) then
     print *, "Error: Could not open output file: ", output_name
     stop
   end if
   write(out_id, *) !skip line
   write(out_id, '(A)') "Discrete Ordinates Method Solution"

   write(out_id, '(A, I0)') "CELL_AVG_PHI  ", size(phi)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i     x_mid            CELL_AVG_PHI"

   x_val = 0
   do ii = 1, Nx
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, (2*x_val+delx(ii))/2, phi(ii) 
       x_val = x_val+ delx(ii)
   end do

   write(out_id, *) 
   write(out_id, '(A, I0)') "PHI_HAT", size(phi_hat)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i     x_mid            PHI_HAT"

   x_val = 0
   do ii = 1, Nx
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, (2*x_val+delx(ii))/2, phi_hat(ii) 
       x_val = x_val+ delx(ii)
   end do


   write(out_id, '(A, I0)') "CELL_AVG_PSI   ", size(psi_avg)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      mu           x_mid            CELL_AVG_PSI"

   do mm = 1, Nq
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       x_val = 0
       do ii = 1, Nx
          write(out_id, '(1I5, 3(1x, ES16.8))') ii, mu(mm), (2*x_val+delx(ii))/2, psi_avg(ii,mm) 
          x_val = x_val+ delx(ii)
       end do
   end do

   write(out_id, '(A, I0)') "CELL_AVG_Jx    ", size(Jx_avg)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      x_mid            CELL_AVG_Jx"

   x_val = 0
   do ii = 1, Nx
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 3(1x, ES16.8))') ii, (2*x_val+delx(ii))/2, Jx_avg(ii) 
       x_val = x_val+ delx(ii)
   end do

   write(out_id, '(A)') "Left face outgoing angular flux"
   write(out_id, '(A)') "     m   mu        outgoiong angular flux"

   do ii = 1, Nbc
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, mu(ii), psi_left(ii) 
   end do

   write(out_id, '(A)') "right face outgoing angular flux"
   write(out_id, '(A)') "     m   mu         outgoiong angular flux"
   do ii = 1, Nbc
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, mu(ii+Nbc), psi_right(ii) 
   end do

   write(out_id, '(A, I0)') "CELL_EDGE_PHI    ", size(phi_edge)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i     x_mid            CELL_EDGE_PHI"

   x_val = 0

   write(out_id, '(1I5, 2(1x, ES16.8))') 1, x_val, phi_edge(1) 

   do ii = 2, Nx+1
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, (x_val+delx(ii-1)), phi_edge(ii) 
       x_val = x_val+ delx(ii-1)
   end do

   write(out_id, '(A, I0)') "CELL_EDGE_PSI    ", size(psi_edge)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      mu           x_mid            CELL_EDGE_PSI"

   do mm = 1, Nq
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       x_val = 0
       write(out_id, '(1I5, 3(1x, ES16.8))') 1, mu(mm), x_val, psi_edge(1,mm) 
       do ii = 2, Nx+1
          write(out_id, '(1I5, 3(1x, ES16.8))') ii, mu(mm), x_val+delx(ii-1), psi_edge(ii,mm) 
          x_val = x_val+ delx(ii-1)
       end do
   end do

   write(out_id, '(A, I0)') "CELL_EDGE_Jx   ", size(Jx_edge)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      x_mid            CELL_EDGE_Jx"

   x_val = 0
   write(out_id, '(1I5, 3(1x, ES16.8))') 1, x_val, Jx_edge(1) 
   do ii = 2, Nx+1
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 3(1x, ES16.8))') ii, x_val+delx(ii-1), Jx_edge(ii) 
       x_val = x_val+ delx(ii-1)
   end do

   write(out_id, '(A, I0)') "SCHEME_RESIDUAL   ", size(scheme_residual)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      x_mid            SCHEME_RESIDUAL"

   x_val = 0

   do ii = 1, Nx
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, (2*x_val+delx(ii))/2, scheme_residual(ii) 
       x_val = x_val+ delx(ii)
   end do

   write(out_id, '(A, I0)') "SCHEME_RESIDUAL_HAT   ", size(scheme_residual_hat)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      x_mid            SCHEME_RESIDUAL_HAT"

   x_val = 0

   do ii = 1, Nx
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, (2*x_val+delx(ii))/2, scheme_residual_hat(ii) 
       x_val = x_val+ delx(ii)
   end do

   write(out_id, '(A, I0)') "ITERATIVE_RESIDUAL   ", size(iterative_residual)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      x_mid            ITERATIVE_RESIDUAL"

   x_val = 0

   do ii = 1, Nx
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, (2*x_val+delx(ii))/2, iterative_residual(ii) 
       x_val = x_val+ delx(ii)
   end do

   write(out_id, '(A, I0)') "ITERATIVE_RESIDUAL_HAT   ", size(iterative_residual_hat)
   ! Printing a descriptive header to match column widths
   write(out_id, '(A)') "     i      x_mid            ITERATIVE_RESIDUAL_HAT"

   x_val = 0

   do ii = 1, Nx
       
       ! Notes on Format:
       ! 3I6: Three integers (i, j, mat), width 6 each
       ! 5(1x, ES16.8): 5 Reals, Scientific, 8 decimals, 1 space padding
       write(out_id, '(1I5, 2(1x, ES16.8))') ii, (2*x_val+delx(ii))/2, iterative_residual_hat(ii) 
       x_val = x_val+ delx(ii)
   end do

   write(out_id, *)  



   write(out_id, *)  
   write(out_id, *) "Stopping criteria, Max inners: ", max_inners 
   write(out_id, "(A, ES16.8)") "Convergence Criterion, e_l < epsilon * (e_l-1/el - 1) with epsilon: ", inner_tol 
   write(out_id, *)  

   if(success_inner) then
     write(out_id,*) "Inner iterations concluded successfully"
   else
     write(out_id, *) "Inner iterations concluded unsuccessfully"
   end if 

   write(out_id, *)  
   write(out_id, "(A, E12.6)") "Number of inners: ", real(inin)
   write(out_id, "(A, ES16.8)") "Achieved convergence for phi: ", error_phi
   write(out_id, "(A, ES16.8)") "Achieved convergence for keff: ", error_keff

   call cpu_time(exec_time)
   exec_time = exec_time - start
   write(out_id, '(A, F16.3, A)') 'Execution Time: ', exec_time, 's'
   write(*, '(A, F16.3, A)') 'Execution Time: ', exec_time, 's'

   close(out_id)

 end subroutine write_phi_to_file 

end module output_mod

