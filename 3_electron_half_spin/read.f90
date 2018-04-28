subroutine read_write

  USE variables 
  implicit none
 
  integer i,j,k,l,m,n,o
  real*8,  dimension(:,:), allocatable  :: Ov, HOv 
  real*8,  dimension(:,:,:,:), allocatable  :: EEOv 
  character st1,st2,st3,st4
  integer nn,ii

  open(unit = 100, file = 'molecule.dat', status = 'old', action = 'read')
      read(100,*) st1, n_el
      read(100,*) st2, sp_deg
      read(100,*) st3, n_bas
      read(100,*) st4, nuc_repul
      allocate ( Ov(n_bas,n_bas),HOv(n_bas,n_bas),EEOv(n_bas,n_bas,n_bas,n_bas),STAT = AllocateStatus)

     close(100)
  open(unit = 101, file = 'ov_integral.dat', status = 'old', action = 'read')

      do i=1,n_bas 
        read(101,*) Ov(i,:)
      end do
    close(101)

  open(unit = 102, file = 'ham_integral.dat', status = 'old', action = 'read') 

      do i=1,n_bas
        read(102,*) HOv(i,:) 
      end do

     close(102)
  open(unit = 103, file = '2e_integral.dat', status = 'old', action = 'read')
      read(103,*) nn
      do ii=1,nn
        read(103,*) i,j,k,l,EEOv(i,k,j,l)
        EEOv(j,k,i,l) = EEOv(i,k,j,l)
        EEOv(i,l,j,k) = EEOv(i,k,j,l)
        EEOv(j,l,i,k) = EEOv(i,k,j,l)
        EEOv(k,i,l,j) = EEOv(i,k,j,l)
        EEOv(k,j,l,i) = EEOv(i,k,j,l)
        EEOv(l,i,k,j) = EEOv(i,k,j,l)
        EEOv(l,j,k,i) = EEOv(i,k,j,l)
      end do      
     close(103)
   call orthogonal(Ov,HOv,EEOv)
!   write(6,*) "Pass Read"
end subroutine read_write
