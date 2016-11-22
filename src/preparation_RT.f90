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
  allocate(dipole_t(0:Nt_iter),norm_t(0:Nt_iter))
  allocate(v_e_new(0:Nx),v_n_new(0:NR))
  allocate(v_e_old(0:Nx,2),v_n_old(0:NR,2))


! Initial condition
  zwfn_e = wfn_e; zwfn_n = wfn_n
  call zwfn_rho
  call mean_field_pot
  v_e_old(:,1) = v_e(:); v_e_old(:,2) = v_e(:)
  v_n_old(:,1) = v_n(:); v_n_old(:,2) = v_n(:)


  write(*,'(A)')'===== Complete preparatin_RT ====================================================='

  return
end subroutine preparation_RT
