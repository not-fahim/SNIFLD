module inner_mod
  use sweep_mod
  use global_vars
  use output_mod
  implicit none
contains

  subroutine inner( inner_tol , & !convergence criteria
                    max_inners, & !maximum number of inners
                    phi,&          !avg scalar flux 
                    phi_old,&      !avg scalar flux old  
                    phi_edge,&     !edge scalar flux
                    phi_hat,&      !scalar flux moment 
                    phi_hat_old,&  !scalar flux moment old
                    Jx_avg,&    ! avg current
                    Jx_edge,&   ! edge current
                    Jx_hat,&    ! spatial moment of current
                    psi_avg,&   ! avg angular flux 
                    psi_edge,&  ! edge angular flux 
                    psi_hat,&   ! spatial moment of angular flux
                    psi_left,&  ! outgoing left face angular flux
                    psi_right,& ! outgoing right angular flux
                    left_bc,&   ! incoming left angular flux
                    right_bc,&  ! incoming right angular flux
                    error_vector,&  !holds difference between scalar flux moment from old iterate
                    Q_i,&          !source
                    Q_hat_i, &       !source spatial moment
                    error_l_list)       !list of L_inf norm for each iteration

    integer, intent(in) :: max_inners      !maximum number of inners, input
    real(kind(0d0)), intent(in) :: inner_tol !convergence criteria, input
  
    real(kind(0d0)), intent(inout) :: phi(:), phi_old(:) 
    real(kind(0d0)), intent(inout) :: phi_edge(:) 
    real(kind(0d0)), intent(inout) :: phi_hat(:), phi_hat_old(:) 

    real(kind(0d0)), intent(inout) :: psi_avg(:,:) 
    real(kind(0d0)), intent(inout) :: psi_edge(:,:) 
    real(kind(0d0)), intent(inout) :: psi_hat(:,:)
    real(kind(0d0)), intent(inout) :: psi_left(:), psi_right(:) 
    real(kind(0d0)), intent(in)    :: left_bc(:), right_bc(:)

    real(kind(0d0)), intent(inout) :: Jx_avg(:) 
    real(kind(0d0)), intent(inout) :: Jx_edge(:) 
    real(kind(0d0)), intent(inout) :: Jx_hat(:)

    !these variable is declared and allocated outside the inner solver 
    !for avoiding allocation in each outer iteration later. 

    real(kind(0d0)), intent(inout) :: error_vector(:) !definded in global, holds difference between the new and old scalar flux

    

    real(kind(0d0)), intent(inout) :: Q_i(:) !holds souce values, 
    real(kind(0d0)), intent(inout) :: Q_hat_i(:) !holds source spatial moments, 

    real(kind(0d0)), intent(inout) :: error_l_list(:) !holds source spatial moments, 

    integer :: i, mat_idx 
    real(kind(0d0)) :: error_l_old, rho_l
    logical                        :: success_inner

    error_l = 1e20 
    error_l_old = 1e20 
    rho_l = 0d0
    

    success_inner = .false. 

    do inin = 1, max_inners 
        phi_old = phi
        phi_hat_old = phi_hat

        !commpute the isotropic souce term for this iteration
        !as the source is isotropc, it is the same for all m
        !so a single value for each cell is used
        do i = 1,Nx
          mat_idx = matmap(i)
          Q_i(i) = q_ext_i(i) + phi_old(i) * (sigma_s(mat_idx) + nu(mat_idx) * sigma_f(mat_idx))
          Q_hat_i(i) = q_ext_hat_i(i) + phi_hat_old(i) * (sigma_s(mat_idx) + nu(mat_idx) * sigma_f(mat_idx))
        end do


        phi = 0d0
        phi_hat = 0d0
        phi_edge = 0d0
        
        Jx_avg = 0d0 
        Jx_edge = 0d0 
        Jx_hat = 0d0
    
        psi_avg = 0d0
        psi_edge = 0d0
        psi_hat = 0d0

        call sweep(   phi,&          !avg scalar flux 
                      phi_edge,&     !edge scalar flux
                      phi_hat,&      !scalar flux moment 

                      Jx_avg,&    ! avg current
                      Jx_edge,&   ! edge current
                      Jx_hat,&    ! spatial moment of current

                      psi_avg,&   ! avg angular flux 
                      psi_edge,&  ! edge angular flux 
                      psi_hat,&   ! spatial moment of angular flux

                      psi_left,&  ! outgoing left face angular flux
                      psi_right,& ! outgoing right angular flux

                      left_bc,&   ! incoming left angular flux
                      right_bc,&  ! incoming right angular flux

                      Q_i,&          !source
                      Q_hat_i)       !source spatial moment


        error_vector = phi - phi_old

        error_l = inf_norm(error_vector)
        rho_l = error_l / error_l_old


        error_l_list(inin) = error_l
        rho_l_list(inin) = rho_l
        if(is_converged(error_l, rho_l, inner_tol, inin)) then
          success_inner = .true. 
          exit  
        end if
        error_l_old = error_l
    end do

    inin = inin-1 

    if(success_inner) then 
        inin = inin+1
        write (*, '(A)')"Complete: inner iterations concluded successfully"
      else
        write(*,'(A)') "Warning: inner iterations concluded unsuccessfully"
    end if

    num_converged_inner = inin 


    write (*, '(A, 1ES16.1)')"number of inner iteration: ", real(inin)
    write(*, *) "Error INF Norm: ", error_l

    !call write_phi_to_file(input_name, start)
    write(*, '(3A)') "Output is saved to: ", trim(input_name), ".out"
  end subroutine inner


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !helper funciton to compute the inf norm
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind(0d0)) function inf_norm(vec)

    real(kind(0d0)), intent(in) :: vec(:) 

    inf_norm = maxval(abs(vec))

  end function inf_norm

  logical function is_converged(error_l, rho_l, epsilon_phi, iter)

    real(kind(0d0)), intent(in) :: error_l, rho_l, epsilon_phi 
    integer, intent(in)         :: iter 
    real(kind(0d0))             :: rhs

    rhs = epsilon_phi * (1/rho_l  - 1)

!   write(*, *) "iteration: ", iter
!   write(*, *) "spectral radius: ", rho_l
!   write(*, *) "rhs for convergence: ", rhs 
!   write(*, *) "lhs for convergence: ", error_l 
!   write(*, *) "" 
!   write(*, *) "" 

    is_converged = .false.
    if(iter<3) return

    if(error_l < rhs) then
        is_converged = .true.
    end if
  end function is_converged
end module inner_mod
