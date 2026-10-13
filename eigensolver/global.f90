module global_vars

  integer                      :: Nx          ! grid dimensions
  real(kind(0d0)), allocatable :: delx(:)     ! grid

  integer                      :: Nq, Nbc     ! quadrature angles, number of BC on each face Nq/2
  real(kind(0d0)), allocatable :: mu(:),  w(:)

  real(kind(0d0)), allocatable :: left_bc(:), right_bc(:) !boundary condition in discrete ordinates

  integer                      :: Nm             ! number of materials
  real(kind(0d0)), allocatable :: sigma_t(:), sigma_s(:), nu(:), sigma_f(:)


  integer, allocatable         :: matmap(:)   ! material map

  ! input source map and calculated source map 
  real(kind(0d0)), allocatable ::  Q_i(:), Q_hat_i(:) !stores the total scattering and ext. src
  real(kind(0d0)), allocatable ::  q_ext_i(:), q_ext_hat_i(:) !stores the ext. src

  !scalar flux map 
  real(kind(0d0)), allocatable ::  phi(:) , phi_old(:)  
  real(kind(0d0)), allocatable ::  phi_hat(:) , phi_hat_old(:), error_vector(:)
  real(kind(0d0)), allocatable ::  phi_edge(:)   

  real(kind(0d0)), allocatable ::  Jx_avg(:) 
  real(kind(0d0)), allocatable ::  Jx_edge(:)   
  real(kind(0d0)), allocatable ::  Jx_hat(:)   

  real(kind(0d0)), allocatable ::  psi_avg(:,:)
  real(kind(0d0)), allocatable ::  psi_hat(:,:)
  real(kind(0d0)), allocatable ::  psi_edge(:,:)   

  real(kind(0d0)), allocatable ::  psi_left(:) , psi_right(:)  ! outgoing angular fluxy
  real(kind(0d0)), allocatable ::  psi_left_old(:) , psi_right_old(:)  ! outgoing angular flux at the boundary
  !scalar flux moment map  

  character(len=8)    :: date_str
  character(len=8)    :: time_str
  character(len=255)  :: input_name 
  real                :: start
  
  !variables for inner iterations
  integer                      :: max_inners 
  integer                      :: inin, success_inner
  real(kind(0d0))              :: inerror, inner_tol
  real(kind(0d0))              :: error_phi, error_phi_old;
  
  !variables for eigen solver
  real(kind(0d0))              :: keff, keff_old;
  real(kind(0d0))              :: fiss_src, fiss_src_old;
  real(kind(0d0))              :: error_keff, error_keff_old;
end module global_vars
