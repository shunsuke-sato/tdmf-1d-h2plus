!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
module global_variables


! parameter
  integer,parameter :: dp = kind(0d0), zp = kind((1d0,1d0))
  complex(zp),parameter :: zI = (0d0,1d0)
  real(dp),parameter :: pi=3.14159265359d0
  
  real(8),parameter :: mass_H = 1836.15267389d0
  real(8),parameter :: mu_e = 2d0*mass_H/(2d0*mass_H+1d0)
  real(8),parameter :: mu_n = 0.5d0*mass_H


! mesh
  integer :: Nx,NR
  real(dp) :: length_x,dx,length_r,dr
  real(dp),allocatable :: xn(:),Rn(:),xRn(:,:)

! GS: wave-function, density, potential and so on
  integer :: Ncg
  real(dp),allocatable :: wfn(:,:), rho(:), v_ext(:), v_int(:,:), v_all(:,:), rho_n(:)
  real(dp),allocatable :: wfn_e(:), wfn_n(:), v_e(:), v_n(:), v0_n(:)
  real(dp),allocatable :: v_ne(:,:), v_en(:,:)
  real(dp),allocatable :: twfn_e(:), twfn_n(:), hwfn_e(:), hwfn_n(:)
  real(dp),allocatable :: twfn_e_b(:), twfn_n_b(:)

! TD
  character(4) :: RT_mode
  integer :: Nt_iter
  real(dp) :: T_calc,dt,kick_mom
  complex(zp),allocatable :: zwfn(:,:)
  real(dp),allocatable :: dipole_t(:),norm_t(:)
  real(dp) :: field_max,field_duration,field_omega
  real(dp) :: field_max_eV_per_AA,field_duration_fs,field_omega_eV
  real(dp),allocatable :: field_t(:)

! temporary
  real(dp),allocatable :: tmp_wfn(:,:),tmp_hwfn(:,:),tmp_wfn_b(:,:)
  complex(zp),allocatable :: ztmp_wfn(:,:),ztmp_hwfn(:,:),ztmp_wfn_b(:,:)

end module global_variables
