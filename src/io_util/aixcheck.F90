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

subroutine AixCheck()

use Fast_IO, only: CtlBlk, FCtlBlk, MxFile, pStat
use Definitions, only: iwp, u6, u5

implicit none
integer(kind=iwp) :: n
logical(kind=iwp) :: is_open
character(len=256) :: filename

!----------------------------------------------------------------------*
! Check if slot in table is available, it should NOT!                  *
!----------------------------------------------------------------------*
do n=1,MxFile
  if (CtlBlk(pStat,n) /= 0) then
    call SysWarnMsg('AixCheck','Active unit: '//FCtlBlk(n),', should have been closed!')
#   ifndef _DEVEL_
    call Abend()
#   endif
  end if
  inquire(unit=n,opened=is_open)
  if (is_open) then
    if ((n /= u6) .and. (n /= u5)) then
      inquire(unit=n,name=filename)
      write(u6,*) 'Fortran file:',n,'(',trim(filename),')  is still open!'
#     ifndef _DEVEL_
      call Abend()
#     endif
    end if
  end if
end do

return

end subroutine AixCheck
