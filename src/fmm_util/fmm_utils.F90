!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module fmm_utils

use fmm_global_paras, only: INTK, REALK, LUPRI, Zero

implicit none
private
! Public procedures
public :: fmm_second, &
          TIMTXT, &
          fmm_quit, &
          fmm_matrix_norm, &
          fmm_distance_between_lines

contains

!-------------------------------------------------------------------------------

function fmm_second()

  implicit none
  real(REALK) :: fmm_second
  integer :: cnt, count_rate, count_max ! note default integer kind for system_clock call
  integer(INTK) :: tmp
  integer(INTK), save :: prt = 0
  integer(INTK), save :: last_count = 0

  ! VV: order changed due to a bug in NAG compiler.
  call system_clock(cnt,count_rate,count_max)

  if (prt == 0) prt = 1  ! first call to fmm_second()
  prt = -prt
  ! must check if cnt has "cycled"
  tmp = cnt
  if (prt == 1) then  ! second of a paired call to fmm_second()
    if (cnt < last_count) tmp = cnt+(count_max-last_count)
  end if
  fmm_second = real(tmp,REALK)/real(count_rate,REALK)
  last_count = cnt

end function fmm_second

!-------------------------------------------------------------------------------

subroutine TIMTXT(TEXTIN,TIMUSD,IUNIT)

  ! TIMTXT based on TIMER by TUH //900709-hjaaj

  implicit none
  character(len=*) :: TEXTIN
  character(len=6) :: AHOUR
  character(len=8) :: ASEC, AMIN
  real(REALK) :: TIMUSD
  integer(INTK) :: ISECND, IUNIT, IHOURS, MINUTE
  character(len=45) :: TEXT

  TEXT = TEXTIN

  ISECND = int(TIMUSD)
  if (ISECND >= 60) then
    MINUTE = ISECND/60
    IHOURS = MINUTE/60
    MINUTE = MINUTE-60*IHOURS
    ISECND = ISECND-3600*IHOURS-60*MINUTE
    if (IHOURS == 1) then
      AHOUR = ' hour '
    else
      AHOUR = ' hours'
    end if
    if (MINUTE == 1) then
      AMIN = ' minute '
    else
      AMIN = ' minutes'
    end if
    if (ISECND == 1) then
      ASEC = ' second '
    else
      ASEC = ' seconds'
    end if
    if (IHOURS > 0) then
      write(IUNIT,100) TEXT,IHOURS,AHOUR,MINUTE,AMIN,ISECND,ASEC
    else
      write(IUNIT,200) TEXT,MINUTE,AMIN,ISECND,ASEC
    end if
  else
    write(IUNIT,300) TEXT,TIMUSD
  end if
  return
  100 format(1x,a,i4,a,i3,a,i3,a)
  200 format(1x,a,i3,a,i3,a)
  300 format(1x,a,f7.2,' seconds')
end subroutine TIMTXT

!-------------------------------------------------------------------------------

subroutine fmm_quit(msg)
  implicit none
  character(len=*) msg
  write(LUPRI,*) msg
  write(LUPRI,*) '>>> FATAL ERROR'
  call Abend()
end subroutine fmm_quit

!------------------------------------------------------------------------------

subroutine fmm_matrix_norm(label,matrix,ndim)

  implicit none
  character(len=*), intent(in) :: label
  integer(INTK), intent(in)    :: ndim
  real(REALK), intent(in)      :: matrix(ndim)

  integer(INTK) :: i
  real(REALK) :: norm

  norm = zero
  do i=1,ndim
    norm = norm+matrix(i)*matrix(i)
  end do
  write(LUPRI,*) 'o fmm_matrix_norm: ',label,' = ',sqrt(norm)

end subroutine fmm_matrix_norm

!------------------------------------------------------------------------------

function fmm_cross_product(a,b)

  implicit none
  real(REALK), intent(in) :: a(3), b(3)
  real(REALK) :: fmm_cross_product(3)

  fmm_cross_product(1) = a(2)*b(3)-a(3)*b(2)
  fmm_cross_product(2) = -(a(1)*b(3)-a(3)*b(1))
  fmm_cross_product(3) = a(1)*b(2)-a(2)*b(1)

end function fmm_cross_product

!-------------------------------------------------------------------------------

function fmm_distance_between_lines(p,q,ep,eq)

  implicit none
  real(REALK), intent(in) :: p(3), q(3), ep(3), eq(3)
  real(REALK) :: fmm_distance_between_lines

  real(REALK) :: d1, d2, n(3)

  n = fmm_cross_product(ep,eq)

  d1 = dot_product(p,n)
  d2 = dot_product(q,n)

  fmm_distance_between_lines = abs(d1-d2)

end function fmm_distance_between_lines

!-------------------------------------------------------------------------------

end module fmm_utils
