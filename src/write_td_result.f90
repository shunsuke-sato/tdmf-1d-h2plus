!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine write_td_results
  use global_variables
  implicit none
  integer :: it


  open(20,file="td_result.out")
  do it = 0,Nt_iter
    write(20,"(999e26.16e3)")it*dt,dipole_t(it),norm_t(it)
  end do
  close(20)

  return
end subroutine write_td_results
