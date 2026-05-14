!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************

subroutine RdTwoEnrg(LU,E2act)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: LU
real(kind=wp), intent(out) :: E2act
integer(kind=iwp) :: istatus
logical(kind=iwp) :: Exists
character(len=80) :: LINE

call OpnFl('INPORB',LU,Exists)
if (.not. Exists) then
  write(u6,*) 'RdTwoEnrg: INPORB not found!'
  call Abend()
end if
rewind(LU)
do
  read(LU,'(A80)',iostat=istatus) Line
  if (istatus /= 0) call Error_quit()
  if (Line(1:22) == '* ACTIVE TWO-EL ENERGY') exit
end do
read(LU,'(ES19.12)',iostat=istatus) E2act
if (istatus /= 0) call Error_quit()

close(LU)

return

contains

subroutine Error_quit()

  call SysWarnFileMsg('RdTwoEnrg','INPORB','Error during reading INPORB\n','Field not there')
  call Abend()

end subroutine Error_quit

end subroutine RdTwoEnrg
