!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine mean_field_pot
  use global_variables
  implicit none
  integer :: ix,iy

  do ix=0,Nx
    v_e(ix) = sum(v_en(:,ix)*rho_n(:))*dr
  end do

  do iy=0,NR
    v_n(iy) = sum(v_ne(:,iy)*rho(:))*dx + v0_n(iy)
  end do

  return
end subroutine mean_field_pot
