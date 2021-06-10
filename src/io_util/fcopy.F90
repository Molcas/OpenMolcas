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

subroutine fcopy(In,Ut,iErr)

character*(*) In, Ut
character*1024 myIn, myUt
integer c_open, c_copy, c_openw, c_close
integer rcIn, rcUt, rc

iErr = 0
if ((len(In) > 1024) .or. (len(Ut) > 1024)) then
  write(6,*) 'Error in fcopy: long filenames'
  iErr = 1
  return
end if
call PrgmTranslate(In,myIn,lIn)
myIn(1+lIn:1+lIn) = char(0)
call PrgmTranslate(Ut,myUt,lUt)
myUt(1+lUt:1+lUt) = char(0)
rcIn = c_open(myIn)
if (rcIn < 0) then
  write(6,*) 'Can not open file ',myIn(1:lIn)
  iErr = 1
  return
end if
rcUt = c_openw(myUt)
if (rcUt < 0) then
  write(6,*) 'Can not open file ',myUt(1:lUt)
  iErr = 1
  return
end if
rc = c_copy(rcIn,rcUt)
if (rc < 0) then
  write(6,*) 'Can not copy file ',myIn(1:lIn)
  iErr = 1
  return
end if
rc = c_close(rcIn)
if (rc < 0) then
  write(6,*) 'Can not close file ',myIn(1:lIn)
  iErr = 1
  return
end if
rc = c_close(rcUt)
if (rc < 0) then
  write(6,*) 'Can not close file ',myUt(1:lUt)
  iErr = 1
  return
end if

return

end subroutine fcopy
