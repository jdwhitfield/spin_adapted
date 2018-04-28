subroutine orthogonal(Ov,HOv,EEOv)
 USE variables
 implicit none
 
 integer i,j,k,l,m,n,o
 integer INFO, LWORK
 real*8  Ov(n_bas,n_bas),HOv(n_bas,n_bas),EEOv(n_bas,n_bas,n_bas,n_bas)  
 real*8  X1(n_bas,n_bas),W(n_bas),WORK(3*n_bas-1),S1(n_bas,n_bas)
 integer a,b,c,d


 allocate ( S(n_bas,n_bas),H(n_bas,n_bas),EE(n_bas,n_bas,n_bas,n_bas), &
          STAT = AllocateStatus)

!-------------------------------------------------------------
! Working in orthogonal basis 

 S1=Ov

 LWORK=3*n_bas-1
 call DSYEV( 'V', 'U', n_bas, Ov, n_bas, W, WORK, LWORK, INFO )

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
!  print*,S(i,:)
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
!write(6,*) "Pass Ortho"
end subroutine 
