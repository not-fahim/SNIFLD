module transport_mod
  use inner_mod 
  use global_vars
  implicit none
contains
  subroutine transport_solver()
    integer :: ii, mm
    allocate (Q_i(Nx), source = 0d0)
    allocate (Q_hat_i(Nx), source = 0d0)

    allocate (psi_left(Nbc), source = 0d0)
    allocate (psi_right(Nbc), source = 0d0)

    allocate (psi_avg(Nx, Nq), source = 0d0)
    allocate (psi_hat(Nx, Nq), source = 0d0)
    allocate (psi_edge(Nx+1, Nq), source = 0d0)


    allocate (phi(Nx), source = 1d0)
    allocate (phi_edge(Nx+1), source = 0d0)
    allocate (phi_hat(Nx), source = 0d0)

    allocate (Jx_avg(Nx), source = 0d0)
    allocate (Jx_edge(Nx+1), source = 0d0)
    allocate (Jx_hat(Nx), source = 0d0)

    allocate (phi_old(Nx), source = 1d0)
    allocate (phi_hat_old(Nx), source = 0d0)
    allocate (error_vector(Nx), source = 0d0)

    keff = 1.00
    keff_old = 1.00
    
    fiss_src = 1.0
    fiss_src_old = 1.0

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
               keff, keff_old,&  !keff to manipulate
               fiss_src, fiss_src_old,&  !keff to manipulate
               Q_i, Q_hat_i)          !source moment variables

  end subroutine transport_solver
end module transport_mod
