module sweep_mod
  use global_vars
  use ldsolver_mod
  implicit none

contains
  subroutine sweep( phi,&          !avg scalar flux 
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

    !variables passed to this subroutine
    real(kind(0d0)), intent(inout) :: phi(:)  
    real(kind(0d0)), intent(inout) :: phi_edge(:) 
    real(kind(0d0)), intent(inout) :: phi_hat(:) 

    real(kind(0d0)), intent(inout) :: psi_avg(:,:) 
    real(kind(0d0)), intent(inout) :: psi_edge(:,:) 
    real(kind(0d0)), intent(inout) :: psi_hat(:,:)
    real(kind(0d0)), intent(inout) :: psi_left(:), psi_right(:) 
    real(kind(0d0)), intent(in) :: left_bc(:), right_bc(:)

    real(kind(0d0)), intent(inout) :: Jx_avg(:) 
    real(kind(0d0)), intent(inout) :: Jx_edge(:) 
    real(kind(0d0)), intent(inout) :: Jx_hat(:)

    real(kind(0d0)), intent(in) :: Q_i(:) 
    real(kind(0d0)), intent(in) :: Q_hat_i(:)
  
    
    !following variables are needed to be passed to ddsolve 
    real(kind(0d0)) :: psi_in_mi, sigma_t_i
    real(kind(0d0)) :: Q_mi, Q_hat_mi
    real(kind(0d0)) :: w_m, mu_m, dx_i  
    real(kind(0d0)) :: psi_mi, psi_hat_mi, psi_out_mi 

    
    integer :: i, mat_idx, m,  count_x

    do m = 1, Nq
    
       mu_m = mu(m)          ! Positive ordinate
       w_m = w(m)
       
       if(mu_m > 0d0) then
          if (left_bc(m) == -1.0d0) then !reflective bc
             psi_in_mi = psi_left(m) 
          else                           !prescribed
             psi_in_mi = left_bc(m)
          end if
    
          i = 1

          !edge values for leftmost
          psi_edge(i, m) = psi_in_mi
          phi_edge(i) = phi_edge(i) + w_m*psi_in_mi
          Jx_edge(i) = Jx_edge(i) + mu_m*w_m*psi_in_mi

          do count_x = 1, Nx
             mat_idx = matmap(i)
             sigma_t_i = sigma_t(mat_idx) 
             Q_mi = Q_i(i)
             Q_hat_mi = Q_hat_i(i)
             dx_i = delx(i) 
    
             call ldsolve(psi_in_mi, Q_mi, Q_hat_mi, sigma_t_i, mu_m, dx_i, &
                          psi_mi, psi_hat_mi, psi_out_mi)

             !avg quantities    
             psi_avg(i,m) = psi_mi 
             phi(i) = phi(i) + w_m * psi_mi    
             Jx_avg(i) = Jx_avg(i) + w_m* mu_m * psi_mi 
             
             !spatial moments
             psi_hat(i,m) = psi_hat_mi
             phi_hat(i) = phi_hat(i) + w_m * psi_hat_mi    
             Jx_hat(i) = Jx_hat(i) + w_m* mu_m * psi_hat_mi 
            
             !edge values
             psi_edge(i+1, m) = psi_out_mi
             phi_edge(i+1) = phi_edge(i+1) + w_m*psi_out_mi
             Jx_edge(i+1) = Jx_edge(i+1) + mu_m*w_m*psi_out_mi

             psi_in_mi = psi_out_mi

             i = i + 1
              
         end do
         psi_right(m) = psi_out_mi
       end if
         
       !store outgoing on right 
    
    
       if(mu_m < 0d0) then
            if (right_bc(m - Nbc) == -1.0d0) then ! we know mu(Nq+n) = -mu(n)
               psi_in_mi = psi_right(m-Nbc)
            else
               psi_in_mi = right_bc(m-Nbc)
            end if
    
            ! Run spatial sweep from i = Nx down to 1
            i = Nx

            !edge values for rightmost
            psi_edge(i+1, m) = psi_in_mi
            phi_edge(i+1) = phi_edge(i+1) + w_m*psi_in_mi
            Jx_edge(i+1) = Jx_edge(i+1) + mu_m*w_m*psi_in_mi

            do count_x = 1, Nx
               mat_idx = matmap(i)
               sigma_t_i = sigma_t(mat_idx) 
               Q_mi = Q_i(i)
               Q_hat_mi = Q_hat_i(i)
               dx_i = delx(i) 
    
               call ldsolve(psi_in_mi, Q_mi, Q_hat_mi, sigma_t_i, mu_m, dx_i, &
                            psi_mi, psi_hat_mi, psi_out_mi)
    
    
                

               !avg quantities    
               psi_avg(i,m) = psi_mi 
               phi(i) = phi(i) + w_m * psi_mi    
               Jx_avg(i) = Jx_avg(i) + w_m* mu_m * psi_mi 
               
               !spatial moments
               psi_hat(i,m) = psi_hat_mi
               phi_hat(i) = phi_hat(i) + w_m * psi_hat_mi    
               Jx_hat(i) = Jx_hat(i) + w_m* mu_m * psi_hat_mi 
            
               !edge values
               psi_edge(i, m) = psi_out_mi
               phi_edge(i) = phi_edge(i) + w_m*psi_out_mi
               Jx_edge(i) = Jx_edge(i) + mu_m*w_m*psi_out_mi

               psi_in_mi = psi_out_mi
               i = i - 1
            end do
            
            psi_left(m-Nbc) = psi_out_mi
       end if
    
    end do
  end subroutine sweep

end module sweep_mod
