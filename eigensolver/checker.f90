module check_mod
    use global_vars
    implicit none

contains

    subroutine check_input()

        !local .......................
        integer :: idx, idm ! a remineder that fortran is case insesative, so they are not i j m
        logical :: flag

        ! delx, dely
        if(Nx ==0 ) stop "Error: Number of cells in x can not be zero"
        if (size(delx) /= Nx) stop "Error: array 'delx' size does not match I."

        do idx = 1, Nx 
            if (delx(idx) <= 0.0d0) stop "Error: delx contains non-positive value."
        end do
        

        ! quadrature checks mu, eta, w
        
        if(Nq==0) stop "Error: number of angles/octant can not be zero"  
        if (size(mu) /= Nq)  stop "Error: Size of 'mu' does not match Nm."
        if (size(w) /= Nq)   stop "Error: Size of 'w' does not match Nq."

        ! Check weights are positive
        do idx = 1, Nq
            if (w(idx) <= 0.0d0) stop "Error: quadrature weight must be positive."
            if (mu(idx) > 1.0d0 ) then
                stop "Error: direction cosines in (mu, eta) must be within [-1, 1]."
            end if

        end do

        !sigma_t, sigma_s

        if(Nm==0) stop "Error: number of materials can not be zero"  
        if (size(sigma_t) /= Nm) stop "Error: size of sigma_t does not match Nm."
        if (size(sigma_s) /= Nm) stop "Error: size of sigma_s does not match Nm."

        do idm = 1, Nm
            if (sigma_t(idm) < 0.0d0) stop "Error: negative total cross-section detected."
            if (sigma_s(idm) < 0.0d0) stop "Error: negative scattering cross-section detected."
            
            if (sigma_s(idm) > sigma_t(idm)) then
                print *, "Error: material", idm, "has Sigma_s > Sigma_t "
                stop
            end if
            if (sigma_t(idm) == 0.0) then
                print *, "Warning: material", idm, "has Sigma_t = 0.0 "
            end if

        end do

        ! boundary check 
        
        if (size(left_bc, 1) /= Nbc) stop "Error: must provide 1 real per Discrete ordinate in left BC"
        if (size(right_bc, 1) /= Nbc) stop "Error: must provide 1 real per Discrete ordinate in right BC"

        flag = .false.
        do idm = 1, Nbc
            if(left_bc(idm) /= -1.0d0) then
              if (flag) stop "Error in LEFT BC: can not be a combination of reflective and other BC type"
              if (left_bc(idm) < 0.0d0) stop "Error in LEFT BC: negative incoming flux is not allowed"

            else 
              flag = .true.
            end if
        end do

        flag = .false.

        do idm = 1, Nbc
            if(right_bc(idm) /= -1.0d0) then
              if (flag) stop "Error in Right BC: can not be a combination of reflective and other BC type"
              if (right_bc(idm) < 0.0d0) stop "Error in RIGHT BC: negative incoming flux is not allowed"

            else 
              flag = .true.
            end if
        end do

        ! matmap, srcmap
        
        if (size(matmap, 1) /= Nx) &
            stop "Error: material map dimensions do not match I."
            
        if (size(q_ext_i, 1) /= Nx ) &
            stop "Error: source map dimensions do not match I."

        if (size(q_ext_hat_i, 1) /= Nx ) &
            stop "Error: source first moment map dimensions do not match I."

        ! Check Nmaterial IDs in map
        do idx = 1, Nx
            if (matmap(idx) < 1 .or. matmap(idx) > Nm) then
                print *, "Matmap Error at cell ", idx 
                stop "Error: material map id is out of range [1, Nm]."
            end if
            
            if (q_ext_i(idx) < 0.0d0) then
                 print *, "srcmap q_ext_i Error at cell ", idx
                 stop "Error: Source term cannot be negative."
            end if
        end do

        print *, "inputs check complete: success"

        if(inner_tol <= 0 ) stop "Error: The inner tolerance must be a positive real"
        if(max_inners <= 0 ) stop "Error: Number of inner iteration must be greater than zero"

    end subroutine check_input

end module check_mod
