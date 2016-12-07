!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dipole(dipole_m,norm_m)
  use global_variables
  implicit none
  real(dp) :: dipole_m, norm_m
  real(dp) :: ss
  integer :: ix,iy

  dipole_m = 0d0; norm_m = 0d0
  do ix = 0,Nx
     dipole_m = dipole_m + xn(ix)*abs(zwfn_e(ix))**2
     norm_m = norm_m + abs(zwfn_e(ix))**2
  end do
  dipole_m = dipole_m*dx
  norm_m = norm_m*dx

  ss = 0d0
  do iy = 0,NR
     ss = ss + abs(zwfn_n(iy))**2
  end do
  ss = ss*dr

  dipole_m = dipole_m * ss
  norm_m = norm_m * ss

  return
end subroutine dipole
