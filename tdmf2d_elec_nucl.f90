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
!=========================================================================================
program main
use global_variables
  implicit none

  write(*,'(A)')'Start qm1d'
  call input
  call mesh

  call preparation_GS
  call Imaginary_time
!  call GS_CG

!  stop
!  call preparation_RT
!  call RT_prop
!
!  call write_results

end program main
!=========================================================================================
subroutine input
  use global_variables
  implicit none

! == input parameter == !                                                                                                                       
  write(*,'(A)')'===== Input parameter ============================================================='
  write(*,'(A)')

  Nx = 400
  length_x = 40d0
  Nr = 400
  length_R = 10d0
! GS
  Ncg = 2400
! TD
  T_calc = 20d0
  dt = 0.01
  Nt_iter = aint(T_calc/dt)+1

  RT_mode='kick' ! kick or gs
  kick_mom = 1e-3

  field_max_eV_per_AA = 1d0
  field_duration_fs = 10d0
  field_omega_eV = 1.55d0

  field_max = field_max_eV_per_AA * (0.529d0/27.2d0)
  field_duration = field_duration_fs/0.02418d0
  field_omega = field_omega_eV/27.2d0

  write(*,'(A,2x,I4)')'Nx =',Nx
  write(*,'(A,2x,e26.16e3)')'length_x =',length_x
  write(*,'(A,2x,e26.16e3)')'T_calc =',T_calc
  write(*,'(A,2x,e26.16e3)')'dt =',dt
  write(*,'(A,2x,I10)')'Nt_iter =',Nt_iter


! temporary array
  allocate(tmp_wfn(0:Nx,0:NR),tmp_hwfn(0:Nx,0:NR),tmp_wfn_b(-4:Nx+4,-4:NR+4))
  allocate(ztmp_wfn(0:Nx,0:NR),ztmp_hwfn(0:Nx,0:NR),ztmp_wfn_b(-4:NR+4,-4:NR+4))
  allocate(twfn_e_b(-4:Nx+4), twfn_n_b(-4:NR+4) )


  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete Input parameter ===================================================' 
  return
end subroutine input
!=========================================================================================
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
!=========================================================================================
subroutine preparation_GS
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp,r,rp,rm

  write(*,'(A)')'===== preparatin_GS =============================================================='
  write(*,'(A)')
  write(*,'(A)')

  allocate(wfn(0:Nx,0:NR), v_ext(0:Nx), v_int(0:Nx,0:NR), v_all(0:Nx,0:NR))
  allocate(rho(0:Nx),rho_n(0:NR))
  allocate(wfn_e(0:Nx), wfn_n(0:NR), v_e(0:Nx), v_n(0:NR), v0_n(0:NR))
  allocate(v_en(0:Nx,0:NR),v_ne(0:NR,0:Nx))
  allocate(twfn_e(0:Nx), twfn_n(0:NR), hwfn_e(0:Nx), hwfn_n(0:NR))




  write(*,'(A)')'=== preparing initial wave-function ===='
! wfn(0,:) == wfn(Nx,:) == 0 
  wfn(:,:) = 0d0

  wfn_e = 0d0;   wfn_n = 0d0
  do ix = 1,Nx-1  
!    call random_number(r); r = r-0.5d0
    wfn_e(ix) = exp(-0.5d0*(xn(ix)-1.0)**2)
  end do
  do ix = 1,NR-1
!    call random_number(r); r = r-0.5d0
     wfn_n(ix) = exp(-0.5d0*(xn(ix)-2.0)**2)  + exp(-0.5d0*(xn(ix)+2.0)**2) 
  end do

! Normalize
  tmp = sum(wfn_e(:)**2)*dx
  wfn_e(:)=wfn_e(:)/sqrt(tmp)
  tmp = sum(wfn_n(:)**2)*dr
  wfn_n(:)=wfn_n(:)/sqrt(tmp)


  write(*,'(A)')'=== preparing external potential ===='
  v_ext = 0d0
!  do ix=0,Nx
!!     v_ext(ix) = 0.5d0*xn(ix)**2
!     v_ext(ix) = -2d0/sqrt(1d0+xn(ix)**2)
!  end do

  write(*,'(A)')'=== preparing interaction potential ===='
  do iy=0,NR
  do ix=0,Nx
     rp = xn(ix) + 0.5d0*Rn(iy)
     rm = xn(ix) - 0.5d0*Rn(iy)
     r = Rn(iy)
!     v_int(ix,iy) = 0d0 !1d0/sqrt(1d0+r**2)
     v_int(ix,iy) = -1d0/sqrt(1d0+rp**2)-1d0/sqrt(1d0+rm**2) +1d0/sqrt(0.03d0+r**2)
     v_en(iy,ix) = -1d0/sqrt(1d0+rp**2)-1d0/sqrt(1d0+rm**2)
     v_ne(ix,iy) = -1d0/sqrt(1d0+rp**2)-1d0/sqrt(1d0+rm**2)
  end do
  end do

  do iy=0,NR
    r = Rn(iy)
    v0_n(iy) = 1d0/sqrt(0.03d0+r**2)
  end do

  write(*,'(A)')'=== preparing total potential ===='

  call wfn_rho
  call mean_field_pot

  do iy=0,NR
  do ix=0,Nx
!     v_all(ix,iy) = v_ext(ix) + v_ext(iy) + v_int(ix,iy)
     v_all(ix,iy) = v_int(ix,iy)
  end do
  end do
  
  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete preparatin_GS ====================================================='

  return
end subroutine preparation_GS
!=========================================================================================
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
!=========================================================================================
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
!=========================================================================================
subroutine GS_CG
  use global_variables
  implicit none
  real(dp),allocatable :: xvec(:,:),pvec(:,:),rvec(:,:)
  real(dp),allocatable :: hxvec(:,:),gvec(:,:),hpvec(:,:)

  real(dp) :: xx,pp,xp,xhx,php,xhp,esp,esp_res,gg,gg0
  real(dp) :: ss,lambda,alpha,beta,aa,bb,cc
  integer :: iorb,iorb_t,ix,iy,iter_cg
  real(dp) :: esp_iter(Ncg),esp_res_iter(Ncg)

  allocate( xvec(0:Nx,0:NR),pvec(0:Nx,0:NR),rvec(0:Nx,0:NR) )
  allocate( hxvec(0:Nx,0:NR),gvec(0:Nx,0:NR),hpvec(0:Nx,0:NR) )

  write(*,'(A)')'===== Ground state calculation ===================================================='
  write(*,'(A)')
  write(*,'(A)')

  xvec(:,:)=wfn(:,:)

  tmp_wfn = xvec; call hpsi; hxvec = tmp_hwfn
  
  xx=sum(xvec(:,:)**2)*dx*dr
  xhx=sum(xvec(:,:)*hxvec(:,:))*dx*dr
  lambda=xhx/xx
  do iter_cg=1,Ncg
     gvec(:,:)=2d0*(hxvec(:,:)-lambda*xvec(:,:))/xx
     
     gg0=sum(gvec(:,:)**2)*dx*dr
     select case(iter_cg)
     case(1)
        pvec(:,:)=-gvec(:,:)
     case default
        beta=gg0/gg
        pvec(:,:)=-gvec(:,:)+beta*pvec(:,:)
     end select
     gg=gg0

     tmp_wfn = pvec; call hpsi; hpvec = tmp_hwfn
        
     pp=sum(pvec(:,:)**2)*dx*dr
     php=sum(pvec(:,:)*hpvec(:,:))*dx*dr
     xp=sum(xvec(:,:)*pvec(:,:))*dx*dr
     xhp=sum(hxvec(:,:)*pvec(:,:))*dx*dr

     aa=php*xp-xhp*pp
     bb=php*xx-xhx*pp
     cc=xhp*xx-xhx*xp
     ss=bb**2-4d0*aa*cc
     if(ss > 0d0)then
        alpha=(-bb+sqrt(ss))/(2d0*aa)
     else
        exit
     end if
     
     xvec(:,:)=xvec(:,:)+alpha*pvec(:,:)

     tmp_wfn = xvec; call hpsi; hxvec = tmp_hwfn
     xx=sum(xvec(:,:)**2)*dx*dr
     xhx=sum(xvec(:,:)*hxvec(:,:))*dx*dr
     lambda=xhx/xx
     esp_iter(iter_cg)=lambda
     esp_res_iter(iter_cg)=sum((hxvec(:,:)-lambda*xvec(:,:))**2)*dx*dr
     
  end do

  xvec(:,:)=xvec(:,:)/sqrt(xx)
  tmp_wfn = xvec; call hpsi; hxvec = tmp_hwfn
  esp=sum(xvec(:,:)*hxvec(:,:))*dx*dr
  esp_res=sum((hxvec(:,:)-esp*xvec(:,:))**2)*dx*dr
  wfn(:,:)=xvec(:,:)
!  if(wfn(Nx/2,Nx/2) < 0d0) wfn(:,:) = -wfn(:,:)


  write(*,'(A)')'esp,     esp_res'
  write(*,'(e16.6e3,3x,e16.6e3)')esp,esp_res


  open(90,file = 'GS_data.out')
  write(90,'(A)')'# iter,  esp(iter),  esp_res(iter)'
  do iter_cg =1,Ncg
     write(90,'(I6,2x,e26.16e3,2x,e26.16e3)')iter_cg,esp_iter(iter_cg),esp_res_iter(iter_cg)
  end do
  close(90)

  open(90,file = 'GS_wfn.out')
  write(90,'(A)')'# x, y, wfn(x,y)'
  do ix =1,Nx
  do iy =1,NR
     write(90,'(100e26.16e3)')xn(ix),Rn(iy),wfn(ix,iy)
  end do
     write(90,*)
  end do
  close(90)

  call wfn_rho
  open(90,file = 'GS_rho_e.out')
  write(90,'(A)')'# x,  rho(x)'
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),rho(ix)
  end do
  close(90)
  write(*,*)'rho',sum(rho(:))*dx


  open(90,file = 'GS_rho_e_R2.0.out')
  write(90,'(A)')'# x,  rho(x)'
  ss = sum(wfn(:,80)**2)*dx
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),wfn(ix,80)**2/ss
  end do
  close(90)

  open(90,file = 'GS_rho_e_R2.5.out')
  write(90,'(A)')'# x,  rho(x)'
  ss = sum(wfn(:,100)**2)*dx
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),wfn(ix,100)**2/ss
  end do
  close(90)

  open(90,file = 'GS_rho_e_R3.0.out')
  write(90,'(A)')'# x,  rho(x)'
  ss = sum(wfn(:,120)**2)*dx
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),wfn(ix,120)**2/ss
  end do
  close(90)


  open(90,file = 'GS_rho_n.out')
  write(90,'(A)')'# x,  rho(x)'
  do ix =0,NR
     write(90,'(100e26.16e3)')Rn(ix),rho_n(ix)
  end do
  close(90)
  write(*,*)'rho',sum(rho_n(:))*dr


  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== End Ground state calculation ================================================'  

  return
end subroutine GS_CG
!=========================================================================================
subroutine wfn_rho
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp

  rho(:)=0d0

!  do ix = 0,Nx
!    tmp = 0d0
!    do iy = 0,Nx
!      tmp = tmp + wfn(ix,iy)**2
!    end do
!    rho(ix) = tmp
!  end do

  rho(:) = wfn_e(:)**2
  rho_n(:) = wfn_n(:)**2


  return
end subroutine wfn_rho
!=========================================================================================
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
!=========================================================================================
