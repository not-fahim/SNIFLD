module transport_mod
  use inner_mod 
  use global_vars
  implicit none
contains
  subroutine transport_solver()
    integer :: ii
    allocate (Q_i(Nx), source = 0d0)
    allocate (Q_hat_i(Nx), source = 0d0)

    allocate (psi_left(Nbc), source = 0d0)
    allocate (psi_right(Nbc), source = 0d0)

    allocate (psi_avg(Nx, Nq), source = 0d0)
    allocate (psi_hat(Nx, Nq), source = 0d0)
    allocate (psi_edge(Nx+1, Nq), source = 0d0)


    allocate (phi(Nx), source = 0d0)
    allocate (phi_edge(Nx+1), source = 0d0)
    allocate (phi_hat(Nx), source = 0d0)

    allocate (Jx_avg(Nx), source = 0d0)
    allocate (Jx_edge(Nx+1), source = 0d0)
    allocate (Jx_hat(Nx), source = 0d0)

    allocate (phi_old(Nx), source = 1d0)
    allocate (phi_hat_old(Nx), source = 0d0)
    allocate (error_vector(Nx), source = 0d0)

    allocate (scheme_residual(Nx), source = 1d0)
    allocate (scheme_residual_hat(Nx), source = 1d0)

    allocate (iterative_residual(Nx), source = 1d0)
    allocate (iterative_residual_hat(Nx), source = 1d0)

    allocate (error_l_list(max_inners), source = -1d0)
    allocate (rho_l_list(max_inners), source = -1d0)

    Q_i = q_ext_i
    Q_hat_i = q_ext_hat_i
    error_vector = phi_old - phi 

    call inner(inner_tol, max_inners,& !iteration paramters  
               phi, phi_old,&          !flux variables
               phi_edge,&
               phi_hat, phi_hat_old,&  !flux moment variables
               Jx_avg, Jx_edge, Jx_hat,&          !current variables
               psi_avg, psi_edge, psi_hat,&          !angular flux variables
               psi_left, psi_right,&          ! outgoing angular flux
               left_bc, right_bc,&          ! incoming angular fluux bcs
               error_vector,&  !flux moment variables
               Q_i, Q_hat_i,&          !source moment variables 
               error_l_list)


    !compute scheme residual
    call compute_residual(phi, phi_old, Jx_edge, scheme_residual)
    call compute_residual_hat(phi_hat, phi_hat_old, Jx_avg, Jx_edge, scheme_residual_hat)

    !compute iterative residual
    call compute_residual(phi, phi, Jx_edge, iterative_residual)
    call compute_residual_hat(phi_hat, phi_hat, Jx_avg, Jx_edge, iterative_residual_hat)


  end subroutine transport_solver
  
  subroutine compute_residual(phi,  phi_old, Jx_edge, residual)
      real(kind(0d0)), intent(in)    :: phi(:), phi_old(:), Jx_edge(:)
      real(kind(0d0)), intent(inout) :: residual(:)
      integer :: ii, mm, mat_idx
    
      do ii = 1, Nx              
          mat_idx = matmap(ii)
          residual(ii) = Jx_edge(ii+1) - Jx_edge(ii) + sigma_t(mat_idx) * phi(ii) * delx(ii) &
                        - (sigma_s(mat_idx) + nu(mat_idx) * sigma_f(mat_idx)) * phi_old(ii) * delx(ii) &
                        - q_ext_i(ii) * delx(ii)
      end do
  end subroutine compute_residual

  subroutine compute_residual_hat(phi_hat,  phi_hat_old, Jx_avg, Jx_edge, residual)
      real(kind(0d0)), intent(in)    :: phi_hat(:), phi_hat_old(:), Jx_edge(:), Jx_avg(:)
      real(kind(0d0)), intent(inout) :: residual(:)

      integer :: ii, mm, mat_idx
    
      do ii = 1, Nx              
          mat_idx = matmap(ii)
          residual(ii) = 3d0 * ( Jx_edge(ii+1) + Jx_edge(ii) - 2d0 * Jx_avg(ii) ) &
                        + sigma_t(mat_idx) * phi_hat(ii) * delx(ii)&
                        - ( sigma_s(mat_idx) + nu(mat_idx) * sigma_f(mat_idx) ) * phi_hat_old(ii) * delx(ii) &
                        - q_ext_hat_i(ii) * delx(ii)
      end do
  end subroutine compute_residual_hat

end module transport_mod
