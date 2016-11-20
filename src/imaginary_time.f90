!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine Imaginary_time
  use global_variables
  implicit none
  integer :: Niter_im = 50000,iter
  real(dp) :: dt_im = 0.005
  real(dp) :: Etot,esp_e,esp_n, res_e,res_n
  real(dp) :: ss
  integer :: ix 


  open(30,file="GS_data.out")
  do iter = 0,Niter_im
    twfn_e = wfn_e; twfn_n = wfn_n
    call hpsi
    esp_e = sum(wfn_e*hwfn_e)*dx; esp_n = sum(wfn_n*hwfn_n)*dr
    res_e = sum( (hwfn_e - esp_e*wfn_e)**2)*dx
    res_n = sum( (hwfn_n - esp_n*wfn_n)**2)*dr
    Etot = esp_e + esp_n - 1.0d0*sum(rho*v_e)*dx
    write(30,"(I7,2x,999e16.6e3)")iter,Etot,esp_e,esp_n,res_e,res_n
    write(*,"(I7,2x,999e16.6e3)")iter,Etot,esp_e,esp_n,res_e,res_n

    wfn_e = wfn_e - dt_im * hwfn_e
    wfn_n = wfn_n - 10d0*dt_im * hwfn_n
    ss = sum(wfn_e**2)*dx; wfn_e = wfn_e/sqrt(ss)
    ss = sum(wfn_n**2)*dr; wfn_n = wfn_n/sqrt(ss)

    call wfn_rho
    call mean_field_pot
    
  end do
  close(30)

  open(90,file = 'GS_wfn_e.out')
  write(90,'(A)')'# x,  rho(x)'
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),wfn_e(ix),rho(ix),v_e(ix)
  end do
  close(90)
  write(*,*)'rho',sum(rho(:))*dx

  open(90,file = 'GS_wfn_n.out')
  write(90,'(A)')'# x,  rho(x)'
  do ix =0,NR
     write(90,'(100e26.16e3)')Rn(ix),wfn_n(ix),rho_n(ix),v_n(ix)
  end do
  close(90)
  write(*,*)'rho_n',sum(rho_n(:))*dr


end subroutine Imaginary_time
