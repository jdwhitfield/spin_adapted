subroutine young_frame_S3
  USE variables
  implicit none

  integer i,j,k,l,m,n,o
  integer, dimension(:,:), allocatable  :: aa,p11,p21
  real*8,  dimension(:), allocatable  :: cof_p11,cof_p21
  integer, dimension(:), allocatable  ::  list_p11,list_p21 
  integer non_deg,deg
  integer n_p11,n_p21

  open(unit = 104, file = 'weyl.dat', status = 'old', action = 'read')
  read(104,*) n_weyl  
  read(104,*) non_deg

   allocate ( aa(n_weyl,n_el),list_p11(n_weyl+1),list_p21(n_weyl+1), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

   n_p11 = non_deg*6+(n_weyl-non_deg)*4
   n_p21 = non_deg*4+(n_weyl-non_deg)*6

   allocate (p11(n_p11,n_el),p21(n_p21,n_el), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
 
   allocate (cof_p11(n_p11),cof_p21(n_p21), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"


   do i = 1, n_weyl
    read(104,*) aa(i,1),aa(i,2),aa(i,3)
   end do
   close(104)
  
!  write(6,*) "Inside Young",non_deg*6+(n_weyl-non_deg)*4,non_deg*4+(n_weyl-non_deg)*6

  list_p11(1) =0
  list_p21(1) =0
!------------------------------------------------------------------
! Generating all the phi_11
do i=1,non_deg
j=((i-1)*6)
j=j+1
cof_p11(j) = 1.d0; p11(j,1)= aa(i,1);p11(j,2)= aa(i,2);p11(j,3)=aa(i,3)
j=j+1
cof_p11(j) =-1.d0; p11(j,1)= aa(i,2);p11(j,2)= aa(i,1);p11(j,3)=aa(i,3)
j=j+1
cof_p11(j) = 0.5d0; p11(j,1)= aa(i,3);p11(j,2)= aa(i,2);p11(j,3)=aa(i,1)
j=j+1
cof_p11(j) = 0.5d0; p11(j,1)= aa(i,1);p11(j,2)= aa(i,3);p11(j,3)=aa(i,2)
j=j+1
cof_p11(j) =-0.5d0; p11(j,1)= aa(i,3);p11(j,2)= aa(i,1);p11(j,3)=aa(i,2)
j=j+1
cof_p11(j) =-0.5d0; p11(j,1)= aa(i,2);p11(j,2)= aa(i,3);p11(j,3)=aa(i,1)
list_p11(i+1) =j
end do

!do i=1,lis_p11(9)
!write(6,*) cof_p11(i),p11(i,:)
!end do

!-----------------------------------------------------------------
! Generating phi_12
do i=1,non_deg
j=((i-1)*4)
j=j+1
cof_p21(j) = sqrt(3.d0)/2.d0; p21(j,1)= aa(i,3);p21(j,2)= aa(i,2);p21(j,3)=aa(i,1)
j=j+1
cof_p21(j) =-sqrt(3.d0)/2.d0; p21(j,1)= aa(i,1);p21(j,2)= aa(i,3);p21(j,3)=aa(i,2)
j=j+1
cof_p21(j) =-sqrt(3.d0)/2.d0; p21(j,1)= aa(i,3);p21(j,2)= aa(i,1);p21(j,3)=aa(i,2)
j=j+1
cof_p21(j) = sqrt(3.d0)/2.d0; p21(j,1)= aa(i,2);p21(j,2)= aa(i,3);p21(j,3)=aa(i,1)
list_p21(i+1) =j
end do


!------------------------------------------------------------------
! Generating all the phi_11
!write(6,*) list_p11(non_deg+1)
j = list_p11(non_deg+1)
do i=non_deg+1,n_weyl
j=j+1
cof_p11(j) =  3.d0/2.d0; p11(j,1)= aa(i,3);p11(j,2)= aa(i,2);p11(j,3)=aa(i,1)
j=j+1
cof_p11(j) = -3.d0/2.d0; p11(j,1)= aa(i,1);p11(j,2)= aa(i,3);p11(j,3)=aa(i,2)
j=j+1
cof_p11(j) =  3.d0/2.d0; p11(j,1)= aa(i,3);p11(j,2)= aa(i,1);p11(j,3)=aa(i,2)
j=j+1
cof_p11(j) = -3.d0/2.d0; p11(j,1)= aa(i,2);p11(j,2)= aa(i,3);p11(j,3)=aa(i,1)
list_p11(i+1) =j
end do

!do i=1,6
!write(6,*) p11(i,:)
!end do

!-----------------------------------------------------------------
! Generating phi_12
!write(6,*) list_p21(non_deg+1)
j = list_p21(non_deg+1)
do i=non_deg+1,n_weyl
j=j+1
cof_p21(j) = sqrt(3.d0); p21(j,1)= aa(i,1);p21(j,2)= aa(i,2);p21(j,3)=aa(i,3)
j=j+1
cof_p21(j) = sqrt(3.d0); p21(j,1)= aa(i,2);p21(j,2)= aa(i,1);p21(j,3)=aa(i,3)
j=j+1
cof_p21(j) =-sqrt(3.d0)/2.d0; p21(j,1)= aa(i,3);p21(j,2)= aa(i,2);p21(j,3)=aa(i,1)
j=j+1
cof_p21(j) =-sqrt(3.d0)/2.d0; p21(j,1)= aa(i,1);p21(j,2)= aa(i,3);p21(j,3)=aa(i,2)
j=j+1
cof_p21(j) =-sqrt(3.d0)/2.d0; p21(j,1)= aa(i,3);p21(j,2)= aa(i,1);p21(j,3)=aa(i,2)
j=j+1
cof_p21(j) =-sqrt(3.d0)/2.d0; p21(j,1)= aa(i,2);p21(j,2)= aa(i,3);p21(j,3)=aa(i,1)
list_p21(i+1) =j
end do

!do i=1,list_p11(n_weyl+1)
! write(6,*) i, cof_p11(i),p11(i,:)
!end do

!do i=1,n_weyl
!write(6,*) list_p11(i),list_p21(i)
!end do

!write(6,*) "Pass Young"
call hamiltonian(cof_p11,cof_p21,p11,p21,list_p11,list_p21,n_p11,n_p21)

end subroutine young_frame_S3
