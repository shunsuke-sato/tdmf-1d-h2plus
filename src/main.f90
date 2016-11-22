!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
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
  call preparation_RT
  call RT_prop
!
!  call write_results

end program main
