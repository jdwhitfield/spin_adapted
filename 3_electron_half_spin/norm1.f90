subroutine norm1(a,p,nm)
USE variables
implicit none
integer a(3),p(3)
real*8 nm
! this subroutine computes the overlap <123|456>
nm=S(a(1),p(1))*S(a(2),p(2))*S(a(3),p(3));
!write(6,*) a,p,S(1,2)
!stop
end subroutine norm1

