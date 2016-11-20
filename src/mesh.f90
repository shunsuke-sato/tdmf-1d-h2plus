!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine mesh
  use global_variables
  implicit none
  integer :: ix, iy 

  write(*,'(A)')'===== Making mesh ================================================================'
  write(*,'(A)')
  write(*,'(A)')

  allocate(xn(0:Nx),Rn(0:NR),xRn(0:Nx,0:NR))
  dx = length_x/dble(Nx)
  dr = length_R/dble(NR)
  
  do ix = 0,Nx
     xn(ix) = dx*dble(ix) - 0.5d0*length_x
  end do

  do ix = 0,NR
    Rn(ix) = dr*dble(ix)
  end do
  
  write(*,'(A)')'===== Complete Making mesh ========================================================'
  return
end subroutine mesh
