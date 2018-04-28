subroutine ham(a,p,v)
USE variables
implicit none
integer a(3),p(3)
real*8  v
! this subroutine computes <123|H|456>

v = H(a(1),p(1))*S(a(2),p(2))*S(a(3),p(3))+ &
    H(a(2),p(2))*S(a(1),p(1))*S(a(3),p(3))+ &
    H(a(3),p(3))*S(a(1),p(1))*S(a(2),p(2))+ &
    EE(a(1),a(2),p(1),p(2))*S(a(3),p(3))+   &
    EE(a(1),a(3),p(1),p(3))*S(a(2),p(2))+   &
    EE(a(2),a(3),p(2),p(3))*S(a(1),p(1));

!    write(6,*) a(3),p(3), H(2,3)
!    stop
end subroutine ham

