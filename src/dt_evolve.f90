!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer,parameter :: Nexp=4
  integer :: iexp 
  real(dp) :: dt_t
  complex(zp) :: zfact

! Use of the (approximated) enforced time-reversal symmetry
! Probably, we need to swithc to the enforced time-reversal symmetry
! employing the predictor-corrector method.

  dt_t = dt
  dt = 0.5d0 * dt


! Propagator = Taylor expantions
  zfact = 1d0
  ztwfn_e = zwfn_e; ztwfn_n = zwfn_n
  do iexp = 1,Nexp
    zfact = zfact*(-zI*dt)/dble(iexp)
    call zhpsi

    zwfn_e = zwfn_e + zfact * zhwfn_e
    zwfn_n = zwfn_n + zfact * zhwfn_n

    ztwfn_e = zhwfn_e; ztwfn_n = zhwfn_n

  end do

  v_e_new(:) = 3d0*v_e(:) -3d0*v_e_old(:,1) + v_e_old(:,2)
  v_n_new(:) = 3d0*v_n(:) -3d0*v_n_old(:,1) + v_n_old(:,2)
  v_e_old(:,2) = v_e_old(:,1); v_e_old(:,1) = v_e(:)
  v_n_old(:,2) = v_n_old(:,1); v_n_old(:,1) = v_n(:)
  v_e = v_e_new; v_n = v_n_new


! Propagator = Taylor expantions
  zfact = 1d0
  ztwfn_e = zwfn_e; ztwfn_n = zwfn_n
  do iexp = 1,Nexp
    zfact = zfact*(-zI*dt)/dble(iexp)
    call zhpsi

    zwfn_e = zwfn_e + zfact * zhwfn_e
    zwfn_n = zwfn_n + zfact * zhwfn_n

    ztwfn_e = zhwfn_e; ztwfn_n = zhwfn_n

  end do

  call zwfn_rho
  call mean_field_pot

  dt = dt_t

end subroutine dt_evolve
