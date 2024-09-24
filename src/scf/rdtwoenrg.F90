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

implicit none
integer LU
real*8 E2act
logical Exist
character(len=80) LINE

call OpnFl('INPORB',LU,Exist)
if (.not. Exist) then
  write(6,*) 'RdTwoEnrg: INPORB not found!'
  call Abend()
end if
rewind(LU)
55 read(LU,'(A80)',end=888,ERR=888) Line
if (Line(1:22) /= '* ACTIVE TWO-EL ENERGY') goto 55
read(LU,'(ES19.12)',err=888,end=888) E2act

close(LU)
return
888 call SysWarnFileMsg('RdTwoEnrg','INPORB','Error during reading INPORB\n','Field not there')
call Abend()

end subroutine RdTwoEnrg
