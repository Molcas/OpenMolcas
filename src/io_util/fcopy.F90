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
interface
  function c_close(FileDescriptor) bind(C,name='c_close_')
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_close
    integer(kind=MOLCAS_C_INT) :: FileDescriptor
  end function c_close
  function c_copy(FileDescriptor1,FileDescriptor2) bind(C,name='c_copy_')
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_copy
    integer(kind=MOLCAS_C_INT) :: FileDescriptor1, FileDescriptor2
  end function c_copy
  function c_open(Path) bind(C,name='c_open_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_open
    character(kind=c_char) :: Path(*)
  end function c_open
  function c_openw(Path) bind(C,name='c_openw_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_openw
    character(kind=c_char) :: Path(*)
  end function c_openw
end interface

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
