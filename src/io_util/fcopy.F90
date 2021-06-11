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

subroutine fcopy(NmIn,NmUt,iErr)

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: NmIn, NmUt
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: lIn, lUt, rc, rcIn, rcUt
character(len=1024) :: myIn, myUt
integer(kind=iwp), external :: c_close, c_copy, c_open, c_openw

iErr = 0
if ((len(NmIn) > 1024) .or. (len(NmUt) > 1024)) then
  write(u6,*) 'Error in fcopy: long filenames'
  iErr = 1
  return
end if
call PrgmTranslate(NmIn,myIn,lIn)
myIn(1+lIn:1+lIn) = char(0)
call PrgmTranslate(NmUt,myUt,lUt)
myUt(1+lUt:1+lUt) = char(0)
rcIn = c_open(myIn)
if (rcIn < 0) then
  write(u6,*) 'Can not open file ',myIn(1:lIn)
  iErr = 1
  return
end if
rcUt = c_openw(myUt)
if (rcUt < 0) then
  write(u6,*) 'Can not open file ',myUt(1:lUt)
  iErr = 1
  return
end if
rc = c_copy(rcIn,rcUt)
if (rc < 0) then
  write(u6,*) 'Can not copy file ',myIn(1:lIn)
  iErr = 1
  return
end if
rc = c_close(rcIn)
if (rc < 0) then
  write(u6,*) 'Can not close file ',myIn(1:lIn)
  iErr = 1
  return
end if
rc = c_close(rcUt)
if (rc < 0) then
  write(u6,*) 'Can not close file ',myUt(1:lUt)
  iErr = 1
  return
end if

return

end subroutine fcopy
