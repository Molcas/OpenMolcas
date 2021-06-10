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

implicit integer(a-z)
#include "ctl.fh"
logical Opened
character name*256

!----------------------------------------------------------------------*
! Check if slot in table is available, it should NOT!                  *
!----------------------------------------------------------------------*
do n=1,MxFile
  if (CtlBlk(pStat,n) /= 0) then
#   ifndef _DEVEL_
    call SysAbendFileMsg('AixCheck',FCtlBlk(n),'Active unit.','Should have been closed!')
#   else
    call SysWarnMsg('AixCheck','Active unit: '//FCtlBlk(n),', should have been closed!')
#   endif
  end if
  inquire(unit=n,Opened=Opened)
  if (Opened) then
    if ((n /= 6) .and. (n /= 5)) then
      inquire(unit=n,name=name)
      write(6,*) 'Fortran file:',n,'(',name(1:index(name,' ')),')  is still open!'
#     ifndef _DEVEL_
      call Abend()
#     endif
    end if
  end if
end do

return

end subroutine AixCheck
