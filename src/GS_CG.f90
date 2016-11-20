!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
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
