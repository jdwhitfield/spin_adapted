module variables
  implicit none  
  integer n_el, sp_deg, n_bas, n_weyl
  real*8  nuc_repul
  real*8,  dimension(:,:), allocatable  :: S, H
  real*8,  dimension(:,:,:,:), allocatable  :: EE 
  integer AllocateStatus
end module
