!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine hpsi
  use global_variables
  implicit none
! finite difference
  real(dp),parameter :: cN0=-205d0/72d0,cN1=8d0/5d0
  real(dp),parameter :: cN2=-1d0/5d0,cN3=8d0/315d0
  real(dp),parameter :: cN4=-1d0/560d0    
  integer :: ix,iy
  real(dp) :: c0,c1,c2,c3,c4
  real(dp) :: c0e,c1e,c2e,c3e,c4e,c0n,c1n,c2n,c3n,c4n
! nine-points formula  
  c0e=-0.5d0*cN0/(dx**2)/mu_e
  c1e=-0.5d0*cN1/(dx**2)/mu_e
  c2e=-0.5d0*cN2/(dx**2)/mu_e
  c3e=-0.5d0*cN3/(dx**2)/mu_e
  c4e=-0.5d0*cN4/(dx**2)/mu_e

  c0n=-0.5d0*cN0/(dr**2)/mu_n
  c1n=-0.5d0*cN1/(dr**2)/mu_n
  c2n=-0.5d0*cN2/(dr**2)/mu_n
  c3n=-0.5d0*cN3/(dr**2)/mu_n
  c4n=-0.5d0*cN4/(dr**2)/mu_n

  twfn_e_b(:)=0d0
  twfn_n_b(:)=0d0
!  tmp_wfn_b(0:Nx,0:NR) = tmp_wfn(0:Nx,0:NR)
  twfn_e_b(0:Nx)= twfn_e(0:Nx)
  twfn_n_b(0:NR)= twfn_n(0:NR)
  twfn_n_b(-1)= twfn_n(1)
  twfn_n_b(-2)= twfn_n(2)
  twfn_n_b(-3)= twfn_n(3)
  twfn_n_b(-4)= twfn_n(4)


  do iy=0,NR
    hwfn_n(iy) = c0n*twfn_n_b(iy) &
      + c1n*( twfn_n_b(iy+1) + twfn_n_b(iy-1) ) &
      + c2n*( twfn_n_b(iy+2) + twfn_n_b(iy-2) ) &
      + c3n*( twfn_n_b(iy+3) + twfn_n_b(iy-3) ) &
      + c4n*( twfn_n_b(iy+4) + twfn_n_b(iy-4) ) 
  end do

  do ix=0,Nx
    hwfn_e(ix) = c0e*twfn_e_b(ix) &
      + c1e*( twfn_e_b(ix+1) + twfn_e_b(ix-1) ) &
      + c2e*( twfn_e_b(ix+2) + twfn_e_b(ix-2) ) &
      + c3e*( twfn_e_b(ix+3) + twfn_e_b(ix-3) ) &
      + c4e*( twfn_e_b(ix+4) + twfn_e_b(ix-4) ) 

  end do

  hwfn_n = hwfn_n + v_n*twfn_n
  hwfn_e = hwfn_e + v_e*twfn_e

  return
end subroutine hpsi
