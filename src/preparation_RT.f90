!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation_RT
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp,r,rp,rm

  write(*,'(A)')'===== preparatin_RT =============================================================='
  write(*,'(A)')
  write(*,'(A)')

  allocate(zwfn(0:Nx,0:NR))
  allocate(ztwfn_e(0:Nx), ztwfn_n(0:NR), zhwfn_e(0:Nx), zhwfn_n(0:NR))

! Initial condition
  zwfn = wfn


  write(*,'(A)')'===== Complete preparatin_RT ====================================================='

  return
end subroutine preparation_RT
