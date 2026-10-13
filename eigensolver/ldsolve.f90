module ldsolver_mod
    implicit none

contains
    subroutine ldsolve( psi_in_mi, &
                        Q_mi, Q_hat_mi, &
                        sigma_t_i, &
                        mu_m,&
                        dx_i,&
                        psi_mi, psi_hat_mi,&
                        psi_out_mi)
    
        !known face averaged incoming angular flux on vertical boundary:
        real(kind(0d0)), intent(in) :: psi_in_mi
    
        !known source within the cell:
        real(kind(0d0)), intent(in) :: Q_mi 

        !known source first spatial moment within the cell:
        real(kind(0d0)), intent(in) :: Q_hat_mi 
        
        !known total cross section within the cell:
        real(kind(0d0)), intent(in) :: sigma_t_i 
        
        !known direction cosines with x 
        real(kind(0d0)), intent(in) :: mu_m

        !known cell dimensions along x 
        real(kind(0d0)), intent(in) :: dx_i 
        
        !unknown cell averaged angular flux:
        real(kind(0d0)), intent(out) :: psi_mi 

        !unknown cell angular flux first moment:
        real(kind(0d0)), intent(out) :: psi_hat_mi 
        
        !unknown face averaged outgoing angular flux on vertical boundary:
        real(kind(0d0)), intent(out) :: psi_out_mi
        
        
        !.............locals..................
        ! the system of eqns looks like:
        ! a1 x + b1 y = c1
        ! a2 x + b2 y = c2
        real(kind(0d0)) :: a1, a2, b1, b2, c1, c2 
 
        if(mu_m > 0) then
          a1 = mu_m + sigma_t_i * dx_i
          b1 = mu_m
          c1 = dx_i * 0.5d0 * Q_mi + psi_in_mi * mu_m

          a2 = -3.0d0 * mu_m
          b2 = 3 * mu_m + sigma_t_i * dx_i
          c2 = 0.5 * dx_i * Q_hat_mi - 3 * mu_m * psi_in_mi
           
          
          psi_mi = (c1*b2 - c2*b1) / (a1*b2 - a2*b1) 
          psi_hat_mi = (a1*c2 - a2*c1) / (a1*b2 - a2*b1) 

          
          psi_out_mi = psi_mi + psi_hat_mi 

        else
          a1 = -mu_m + sigma_t_i * dx_i
          b1 = mu_m
          c1 = dx_i * 0.5d0 * Q_mi - psi_in_mi * mu_m

          a2 = -3.0d0 * mu_m
          b2 = -3 * mu_m + sigma_t_i * dx_i
          c2 = 0.5 * dx_i * Q_hat_mi - 3 * mu_m * psi_in_mi
           
          
          psi_mi = (c1*b2 - c2*b1) / (a1*b2 - a2*b1) 
          psi_hat_mi = (a1*c2 - a2*c1) / (a1*b2 - a2*b1) 

          
          psi_out_mi = psi_mi - psi_hat_mi 

        endif
    
    end subroutine ldsolve
end module ldsolver_mod
