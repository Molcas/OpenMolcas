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
! Copyright (C) 1990,1991, Per-Olof Widmark                            *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                             A I X - I / O                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixOpn(Handle,FileName)                                           *
!                                                                      *
! A file is opened for read/write operations. If the file does not     *
! exist it is automatically created.                                   *
!                                                                      *
! Input:  FileName - A character string specifying a complete path     *
!                    name or relative path name for the file. The      *
!                    name must be shorter than 128 characters.         *
!                                                                      *
! Output: Handle   - When a file is sucessfully opened, a unique file  *
!                    id is supplied by the routine. This is used for   *
!                    specifying the file to other routines.            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          S&TC, ACIS, IBM Sweden                                      *
! Written: November 1990                                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
! 911021 - If return code is 13 when opening file, try to open it as a *
!          read only file. Per-Olof Widmark.                           *
! 931217 - Flags obtained from system calls.                           *
! 010710 - Change using iRand to incremental number                    *
! 120304 - Use handle for the passing of 'Lu' value needed by FiM      *
!                                                                      *
!***********************************************************************

function AixOpn(handle,filename,translate)

use Fast_IO, only: CtlBlk, eBlNme, eTlFn, eTmF, FCtlBlk, MxFile, pDesc, pHndle, pStat, pWhere
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: AixOpn
integer(kind=iwp), intent(inout) :: handle
character(len=*), intent(in) :: filename
logical(kind=iwp), intent(in) :: translate
integer(kind=iwp) :: desc, ltmp, n, nFile, NVV = 666, rc
character(len=80) :: ErrTxt
character(len=256) :: tmp, tmp1
integer(kind=iwp), external :: StrnLn
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
  function c_open(Path) bind(C,name='c_open_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_open
    character(kind=c_char) :: Path(*)
  end function c_open
end interface

!----------------------------------------------------------------------*
! Entry to AixOpn                                                      *
!----------------------------------------------------------------------*
AixOpn = 0
!----------------------------------------------------------------------*
! Check if slot in table is available                                  *
!----------------------------------------------------------------------*
n = 1
do
  if (CtlBlk(pStat,n) == 0) exit
  n = n+1
  if (n > MxFile) then
    AixOpn = eTmF
    call SysWarnMsg('Aixopn','Too many opened files\n','try to increase MxFile')
    return
  end if
end do
nFile = n
!----------------------------------------------------------------------*
! Strip file name and append string terminator                         *
!----------------------------------------------------------------------*
n = len(filename)
do
  if (filename(n:n) /= ' ') exit
  n = n-1
  if (n <= 0) then
    AixOpn = eBlNme
    return
  end if
end do
n = n+1
if (n >= len(tmp)) then
  AixOpn = eTlFn
  return
end if
tmp = filename
tmp(n:n) = char(0)
!----------------------------------------------------------------------*
! Attempt to open file.                                                *
!----------------------------------------------------------------------*
rc = 0
tmp1 = tmp
ltmp = StrnLn(tmp1)
if (translate) then
  call PrgmTranslate(tmp1,tmp,ltmp)
end if
!write(u6,*) 'DEBUG AIXOPN: what would happen ',tmp(1:ltmp),'<'
tmp = tmp(1:ltmp)

!write(u6,*) 'in=',tmp1
!write(u6,*) 'len=',ltmp
!write(u6,*) 'res=',tmp
tmp(ltmp+1:ltmp+1) = char(0)
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
if (handle == 0) then
  rc = c_open(tmp)
else
  rc = c_FimOpen(tmp,handle)
  if (rc < 0) then
    rc = -rc
    AixOpn = eFiMFo
  end if
end if
if (handle == 0) then
#else
  rc = c_open(tmp)
#endif
  if (rc < 0) then
    rc = AixErr(ErrTxt)
    call SysWarnFileMsg('AixOpn',filename,'MSG: open',ErrTxt)
    call SysPutsEnd()
    call Abend()
  end if
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
end if
#endif
desc = rc
!----------------------------------------------------------------------*
! Attempt sucessful, update control blocks.                            *
!----------------------------------------------------------------------*
!handle = iRand()
NVV = NVV+100
handle = NVV
CtlBlk(pHndle,nFile) = handle
CtlBlk(pDesc,nFile) = desc
CtlBlk(pStat,nFile) = 1
CtlBlk(pWhere,nFile) = 0
FCtlBlk(nFile) = filename
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixOpn
