program half_spin
implicit none
integer i,j,k,l,a,b,c,d,ii,nn
integer n_weyl,n_bas
integer INFO, LWORK,LWORK_f
real*8  val
real*8  Xsymm(3,3)
real*8  E_nuc
integer AllocateStatus

real*8,  dimension(:), allocatable  :: W_f,WORK_f
real*8,  dimension(:), allocatable  :: W,WORK,e1
real*8,  dimension(:,:), allocatable  :: H_f,Ov, HOv,S,H,S1,X1
real*8,  dimension(:,:,:,:), allocatable  :: EEOv,EE
integer, dimension(:,:), allocatable  :: aa,p11,p12


open(unit = 100, file = 'weyl.dat', status = 'old', action = 'read')
 read(100,*) n_weyl

allocate ( aa(n_weyl,3), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

allocate ( W_f(n_weyl),WORK_f(3*n_weyl-1),H_f(n_weyl,n_weyl), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

allocate (p11(n_weyl*6,4),p12(n_weyl*4,4), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

 do i = 1, n_weyl
  read(100,*) aa(i,1),aa(i,2),aa(i,3)
 end do
 close(100)

! do i = 1,8
! write(6,*) a(i,1),a(i,2),a(i,3)
! end do
!

!------------------------------------------------------------------
! Generating all the phi_11
do i=1,n_weyl
j=((i-1)*6)
j=j+1
p11(j,1) = 2; p11(j,2)= aa(i,1);p11(j,3)= aa(i,2);p11(j,4)=aa(i,3)
j=j+1
p11(j,1) =-2; p11(j,2)= aa(i,2);p11(j,3)= aa(i,1);p11(j,4)=aa(i,3)
j=j+1
p11(j,1) = 1; p11(j,2)= aa(i,3);p11(j,3)= aa(i,2);p11(j,4)=aa(i,1)
j=j+1
p11(j,1) = 1; p11(j,2)= aa(i,1);p11(j,3)= aa(i,3);p11(j,4)=aa(i,2)
j=j+1
p11(j,1) = -1; p11(j,2)= aa(i,3);p11(j,3)= aa(i,1);p11(j,4)=aa(i,2)
j=j+1
p11(j,1) = -1;p11(j,2)= aa(i,2);p11(j,3)= aa(i,3);p11(j,4)=aa(i,1)
end do

!do i=1,6
!write(6,*) p11(i,:)
!end do

!-----------------------------------------------------------------
! Generating phi_12
do i=1,n_weyl
j=((i-1)*4)
j=j+1
p12(j,1) =1; p12(j,2)= aa(i,3);p12(j,3)= aa(i,2);p12(j,4)=aa(i,1)
j=j+1
p12(j,1) =1; p12(j,2)= aa(i,1);p12(j,3)= aa(i,3);p12(j,4)=aa(i,2)
j=j+1
p12(j,1) =-1; p12(j,2)= aa(i,3);p12(j,3)= aa(i,1);p12(j,4)=aa(i,2)
j=j+1
p12(j,1) =-1; p12(j,2)= aa(i,2);p12(j,3)= aa(i,3);p12(j,4)=aa(i,1)

end do

!do i=1,4
!write(6,*) p12(i,:)
!end do

!----------------------------------------------------------------------
! loading the integrals - Overlap , Core Hamiltonian, Nuclear repulsion
!                         two electron integral (input is expected to be in chemist notion)

open(unit = 101, file = 'ham_ov.dat', status = 'old', action = 'read')
 read(101,*) n_bas
 allocate ( Ov(n_bas,n_bas),HOv(n_bas,n_bas),EEOv(n_bas,n_bas,n_bas,n_bas), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
 allocate ( S(n_bas,n_bas),H(n_bas,n_bas),EE(n_bas,n_bas,n_bas,n_bas),S1(n_bas,n_bas), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***" 


 do i = 1,n_bas
  read(101,*) Ov(i,:)
!  write(6,*) Ov(i,:)
 end do

 do i = 1,n_bas
  read(101,*) HOv(i,:)
!  write(6,*) HOv(i,:)
 end do
 read(101,*) E_nuc
 close(101)


open(unit = 102, file = '2_ele.dat', status = 'old', action = 'read')

 allocate ( e1(n_bas**4), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  read(102,*) nn
do ii=1,nn
 read(102,*) i,j,k,l,EEOv(i,k,j,l)
 EEOv(j,k,i,l) = EEOv(i,k,j,l)
 EEOv(i,l,j,k) = EEOv(i,k,j,l)
 EEOv(j,l,i,k) = EEOv(i,k,j,l)
 EEOv(k,i,l,j) = EEOv(i,k,j,l)
 EEOv(k,j,l,i) = EEOv(i,k,j,l)
 EEOv(l,i,k,j) = EEOv(i,k,j,l)
 EEOv(l,j,k,i) = EEOv(i,k,j,l)
end do

S1=Ov

 allocate ( W(n_bas),WORK(3*n_bas-1), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
LWORK=3*n_bas-1
call DSYEV( 'V', 'U', n_bas, Ov, n_bas, W, WORK, LWORK, INFO )

! do i=1,3
!   write(6,*)S1(i,:)
! end do

!-------------------------------------------------------------
! Working in orthogonal basis 

allocate ( X1(n_bas,n_bas), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
 do i=1,n_bas
    do  j=1,n_bas
      X1(i,j)=Ov(i,j)/sqrt(W(j))
    end do 
 end do 

! Symmetric Orthogonalization
! do i=1,3
!    do  j=1,3
!      X1(i,j)=Ov(i,j)/sqrt(W(j))
!    end do 
! end do 
!Xsymm=matmul(X1,transpose(Ov))
!
!X1=Xsymm
!X1=0.d0

!do i=1,3
!X1(i,i)=1.d0
!enddo 


S=matmul(transpose(X1),matmul(S1,X1))
H=matmul(transpose(X1),matmul(HOv,X1))

!do i =1,n_bas
! print*,S(i,:)
!end do

do i=1,n_bas
  do j=1,n_bas
    do k=1,n_bas
      do l=1,n_bas

          EE(i,j,k,l)=0.d0

         do a=1,n_bas
          do b=1,n_bas
           do c=1,n_bas
            do d=1,n_bas
              EE(i,j,k,l)=EE(i,j,k,l)+X1(a,i)*X1(b,j)*X1(c,k)*X1(d,l)*EEOv(a,b,c,d)
            end do 
           end do 
          end do 
         end do 

       end do 
     end do 
   end do 
 end do 


!--------------------------------------------------------------------------------
! computing the Hessian matrix

do i=1,n_weyl
    do j=1,n_weyl
       call one_element(p11(6*(i-1)+1:6*(i-1)+6,:),p12(4*(i-1)+1:4*(i-1)+4,:), &
                        p11(6*(j-1)+1:6*(j-1)+6,:),p12(4*(j-1)+1:4*(j-1)+4,:),val,S,H,EE,n_bas)

      H_f(i,j)=val
      
      if(abs(H_f(i,j)).lt.1.D-8) then
       H_f(i,j)=0.d0
      end if
!      write(6,*) H_f(i,j)
    end do
end do

do i=1,n_weyl
! write(23,10)H_f(i,:)
enddo 

!--------------------------------------------------------------------------------
!Digonalizing the matrix

LWORK_f=3*n_weyl-1 
call DSYEV( 'V', 'U', n_weyl, H_f, n_weyl, W_f, WORK_f, LWORK_f, INFO )

 write(6,*) "Spectrum for Sz=1/2"
do i=1,n_weyl
 write(6,10)W_f(i)+E_nuc
10 format(8f13.8)
enddo


end program half_spin


subroutine one_element(a1,b1,a2,b2,val,S,H,EE,n_bas)
integer i,j,k,l,n_bas
integer a1(6,4),a2(6,4),b1(4,4),b2(4,4)
integer aa1(6,3),aa2(6,3),bb1(4,3),bb2(4,3)
real*8 v,val,val1,val2,nm,norm_aa1,norm_aa2,norm_bb1,norm_bb2
real*8  S(n_bas,n_bas),H(n_bas,n_bas),EE(n_bas,n_bas,n_bas,n_bas)

!--------------------------------------------------------------------------------
! This subroutines compute 
!  part <PHI(i)|H|PHI(j)>  where 
! PHI(i)=phi_11(i)*alpha_1 + phi_12(i)*alpha_2
! PHI(i)=phi_11(j)*alpha_1 + phi_12(j)*alpha_2
! and <alpha_1|alpha_2> = 0 i.e orthogonal
!--------------------------------------------------------------------------------

aa1=a1(:,2:4)   ! storing the components of phi_11(i) 
aa2=a2(:,2:4)   ! storing the components of phi_12(j)


 norm_aa1=0.d0;
 norm_aa2=0.d0;

! Computing norm square of phi_11(i)

 do i=1,6
  do k=1,6
  call norm1(aa1(i,:),aa1(k,:),nm,S,H,EE,n_bas)
  norm_aa1=norm_aa1+a1(i,1)*a1(k,1)*nm/4.d0
  end do
 end do

! Computing norm square of phi_11(j)

 do i=1,6
  do k=1,6
   call norm1(aa2(i,:),aa2(k,:),nm,S,H,EE,n_bas)
   norm_aa2=norm_aa2+a2(i,1)*a2(k,1)*nm/4.d0
  end do
 end do

! Computing the <phi_11(i)|H|phi_11(j)>

val1=0
 do i=1,6
   do j=1,6
       call ham(aa1(i,:),aa2(j,:),v,S,H,EE,n_bas)
       val1=val1+a1(i,1)*a2(j,1)*v/4.d0
   end do
 end do



bb1=b1(:,2:4);    ! storing the components of phi_12(i)
bb2=b2(:,2:4);    ! storing the components of phi_12(j)


 norm_bb1=0.d0;
 norm_bb2=0.d0;

! Computing norm square of phi_12(i)

 do i=1,4
  do k=1,4
   call norm1(bb1(i,:),bb1(k,:),nm,S,H,EE,n_bas)
   norm_bb1=norm_bb1+b1(i,1)*b1(k,1)*nm*3.d0/4.d0
  end do 
 end do 

! Computing norm square of phi_12(j)

 do i=1,4
  do k=1,4
  call norm1(bb2(i,:),bb2(k,:),nm,S,H,EE,n_bas)
  norm_bb2=norm_bb2+b2(i,1)*b2(k,1)*nm*3.d0/4.d0
  end do
 end do

! Computing <phi_12(i)|H|phi_12(j)>

val2=0;
  do i=1,4
   do j=1,4
       call ham(bb1(i,:),bb2(j,:),v,S,H,EE,n_bas)
       val2=val2+b1(i,1)*b2(j,1)*v*3.d0/4.d0
   end do
 end do


! This is more crucial step 
!
!1. This assuming the alpha1 and alpha2 are already orthonormal and normalizing it with norms
!   of individual part
!
   val= (val1+val2)/(sqrt((norm_aa1+norm_bb1))*sqrt((norm_aa2+norm_bb2)))

!2. This assuming the alpha1 and alpha2 are already orthonormal and using the direct expression
!   the book by Ruben Paunz
!
! val=(val1+val2)/(6.d0)

!3. This assuming the alpha1 and alpha2 are not normalized and get constructed from conjugate young
!   tabulae
!   alpha1= sqrt(3)[beta*alpha*alpha - alpha*alpha*beta]
!   alpha2= [2*alpha*alpha*beta - beta*lpha*alpha - alpha*beta*alpha ]
!   <alpha1|alpha1> = 6 ; <alpha2|alpha2>=6
!
!  val= (val1*6.d0+val2*3.d0)/(sqrt((norm_aa1*6.d0+norm_bb1*3.d0))*sqrt((norm_aa2*6.d0+norm_bb2*3.d0)))
end subroutine one_element 


subroutine ham(a,p,v,S,H,EE,n_bas)
integer n_bas
real*8  S(n_bas,n_bas),H(n_bas,n_bas),EE(n_bas,n_bas,n_bas,n_bas)
integer a(3),p(3)
real*8  v
! this subroutine computes <123|H|456>

v = H(a(1),p(1))*S(a(2),p(2))*S(a(3),p(3))+ &
    H(a(2),p(2))*S(a(1),p(1))*S(a(3),p(3))+ &
    H(a(3),p(3))*S(a(1),p(1))*S(a(2),p(2))+ &
    EE(a(1),a(2),p(1),p(2))*S(a(3),p(3))+   &
    EE(a(1),a(3),p(1),p(3))*S(a(2),p(2))+   &
    EE(a(2),a(3),p(2),p(3))*S(a(1),p(1));

end subroutine ham

subroutine norm1(a,p,nm,S,H,EE,n_bas)

integer a(3),p(3),n_bas
real*8  S(n_bas,n_bas),H(n_bas,n_bas),EE(n_bas,n_bas,n_bas,n_bas)
real*8 nm
! this subroutine computes the overlap <123|456>
nm=S(a(1),p(1))*S(a(2),p(2))*S(a(3),p(3));

end subroutine norm1


