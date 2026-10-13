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
                    keff, keff_old,&  !keff to manipulate
                    fiss_src, fiss_src_old,&  !keff to manipulate
                    Q_i,&          !source
                    Q_hat_i)       !source spatial moment

    integer, intent(in) :: max_inners      !maximum number of inners, input
    real(kind(0d0)), intent(in) :: inner_tol !convergence criteria, input
  
    real(kind(0d0)), intent(inout) :: phi(:), phi_old(:) 
    real(kind(0d0)), intent(inout) :: phi_edge(:) 
    real(kind(0d0)), intent(inout) :: phi_hat(:), phi_hat_old(:) 

    real(kind(0d0)), intent(inout) :: psi_avg(:,:) 
    real(kind(0d0)), intent(inout) :: psi_edge(:,:) 
    real(kind(0d0)), intent(inout) :: psi_hat(:,:)
    real(kind(0d0)), intent(inout) :: psi_left(:), psi_right(:) 
    real(kind(0d0)), intent(in) :: left_bc(:), right_bc(:)

    real(kind(0d0)), intent(inout) :: Jx_avg(:) 
    real(kind(0d0)), intent(inout) :: Jx_edge(:) 
    real(kind(0d0)), intent(inout) :: Jx_hat(:)

    real(kind(0d0)), intent(inout) :: keff, keff_old
    real(kind(0d0)), intent(inout) :: fiss_src, fiss_src_old

    !these variable is declared and allocated outside the inner solver 
    !for avoiding allocation in each outer iteration later. 

    real(kind(0d0)), intent(inout) :: error_vector(:) !definded in global, holds difference between the new and old scalar flux

    

    real(kind(0d0)), intent(inout) :: Q_i(:) !holds souce values, 
    real(kind(0d0)), intent(inout) :: Q_hat_i(:) !holds source spatial moments, 

    integer :: i, mat_idx 
    real(kind(0d0)) :: abs_error
    

    success_inner = 0

    do inin = 1, max_inners 

        keff_old = keff
        fiss_src_old = fiss_src
        phi_old = phi
        phi_hat_old = phi_hat

        !commpute the isotropic souce term for this iteration
        !as the source is isotropc, it is the same for all m
        !so a single value for each cell is used
        do i = 1,Nx
          mat_idx = matmap(i)
          Q_i(i) = q_ext_i(i) + &
                   phi_old(i) * (sigma_s(mat_idx) + nu(mat_idx) * sigma_f(mat_idx) / keff_old)
          Q_hat_i(i) = q_ext_hat_i(i) + &
                      phi_hat_old(i) * (sigma_s(mat_idx) + nu(mat_idx) * sigma_f(mat_idx)/ keff_old)
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

        fiss_src = compute_fission_source(phi)

        keff = keff_old * fiss_src / fiss_src_old  
        error_keff = keff - keff_old

        error_vector = phi - phi_old
        error_phi= inf_norm(error_vector)

        if( is_converged(error_phi, error_phi_old, inner_tol, inin) &
            .and. is_converged(error_keff, error_keff_old, inner_tol, inin)  ) then
          success_inner = 1
          exit  
        end if

        error_phi_old = error_phi
        error_keff_old = error_keff
    end do

    inin = inin-1 

    if(success_inner == 1) then 
        inin = inin+1
        write (*, '(A)')"Complete: inner iterations concluded successfully"
      else
        write(*,'(A)') "Warning: inner iterations concluded unsuccessfully"
    end if


    write (*, '(A, 1ES16.1)')"number of inner iteration: ", real(inin)
    write(*, *) "Error INF Norm Phi: ", error_phi
    write(*, *) "Error in Keff: ", error_keff
    write(*, *) "Result Keff: ", keff

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !helper function to compute fission rate from scalar flx
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind(0d0)) function compute_fission_source(phi)

    real(kind(0d0)), intent(in) :: phi(:) 
    integer                     :: ii, mat_idx

    compute_fission_source = 0.0

    do ii = 1, Nx
        mat_idx = matmap(ii)
        compute_fission_source = compute_fission_source + phi(ii) * nu(mat_idx) * sigma_f(mat_idx)
    end do

  end function compute_fission_source
    

  logical function is_converged(error_l, error_l_old, epsilon_phi, iter)

    real(kind(0d0)), intent(in) :: error_l, error_l_old, epsilon_phi 
    real(kind(0d0))             :: rho_l, rhs
    integer                     :: iter 
    


    rho_l = error_l / error_l_old
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
