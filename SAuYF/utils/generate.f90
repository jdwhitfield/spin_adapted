program generate
implicit none
integer i,j,k,l,bas,ii
real*8 a,b,c,d,intg(3,3,3,3)

!write(6,*) "Number of contracted basis"
!read(*,*) bas

!intg = 0.d0

do ii=1,21
 read (*,*) i,j,k,l,intg(i,j,k,l) 
 intg(i,j,l,k) = intg(i,j,k,l)
 intg(j,i,k,l) = intg(i,j,k,l)
 intg(j,i,l,k) = intg(i,j,k,l)
 intg(k,l,i,j) = intg(i,j,k,l)
 intg(l,k,i,j) = intg(i,j,k,l)
 intg(k,l,j,i) = intg(i,j,k,l)
 intg(l,k,j,i) = intg(i,j,k,l)
end do

open (unit = 22, file = "2_ele_chemist.dat")
open (unit = 23, file = "2_ele.dat")
open (unit = 24, file = "2_ele_physicist.dat")

do i=1,3
 do j=1,3
  do k=1,3
   do l=1,3
    write(22,*) i,j,k,l,intg(i,j,k,l)
    write(23,10) intg(i,k,j,l)
    write(24,*) i,k,j,l, intg(i,k,j,l)
    10 format(8f13.11)
   end do
  end do
 end do
end do
close(22)
close(23)
end program generate
