!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine zhpsi
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

  ztwfn_e_b(:)=0d0
  ztwfn_n_b(:)=0d0
!  tmp_wfn_b(0:Nx,0:NR) = tmp_wfn(0:Nx,0:NR)
  ztwfn_e_b(0:Nx)= ztwfn_e(0:Nx)
  ztwfn_n_b(0:NR)= ztwfn_n(0:NR)
  ztwfn_n_b(-1)= ztwfn_n(1)
  ztwfn_n_b(-2)= ztwfn_n(2)
  ztwfn_n_b(-3)= ztwfn_n(3)
  ztwfn_n_b(-4)= ztwfn_n(4)


  do iy=0,NR
     zhwfn_n(iy) = c0n*ztwfn_n_b(iy) &
      + c1n*( ztwfn_n_b(iy+1) + ztwfn_n_b(iy-1) ) &
      + c2n*( ztwfn_n_b(iy+2) + ztwfn_n_b(iy-2) ) &
      + c3n*( ztwfn_n_b(iy+3) + ztwfn_n_b(iy-3) ) &
      + c4n*( ztwfn_n_b(iy+4) + ztwfn_n_b(iy-4) ) 
  end do

  do ix=0,Nx
     zhwfn_e(ix) = c0e*ztwfn_e_b(ix) &
      + c1e*( ztwfn_e_b(ix+1) + ztwfn_e_b(ix-1) ) &
      + c2e*( ztwfn_e_b(ix+2) + ztwfn_e_b(ix-2) ) &
      + c3e*( ztwfn_e_b(ix+3) + ztwfn_e_b(ix-3) ) &
      + c4e*( ztwfn_e_b(ix+4) + ztwfn_e_b(ix-4) ) 

  end do

  zhwfn_n = zhwfn_n + v_n*ztwfn_n
  zhwfn_e = zhwfn_e + v_e*ztwfn_e

  return
end subroutine zhpsi
