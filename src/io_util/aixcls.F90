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
! Copyright (C) 1990, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                             A I X - I / O                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixCls(Handle)                                                    *
!                                                                      *
! A file is closed.                                                    *
!                                                                      *
! Input:  Handle   - This is the unique file identifier associated     *
!                    with the file. It is created by AixOpn, and must  *
!                    be used on subsequent references to the file.     *
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
!                                                                      *
!***********************************************************************

function AixCls(handle)

use Fast_IO, only: CtlBlk, eNtOpn, FCtlBlk, MxFile, pDesc, pHndle, pStat, pWhere
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: AixCls
integer(kind=iwp), intent(in) :: handle
integer(kind=iwp) :: desc, n, nFile, rc
character(len=80) :: ErrTxt
integer(kind=iwp), external :: c_close
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
end interface

!----------------------------------------------------------------------*
! Entry to AixCls                                                      *
!----------------------------------------------------------------------*
AixCls = 0
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
do
  if (CtlBlk(pHndle,n) == handle) exit
  n = n+1
  if (n > MxFile) then
    AixCls = eNtOpn
    return
  end if
end do
nFile = n
desc = CtlBlk(pDesc,nFile)
!----------------------------------------------------------------------*
! Close file                                                           *
!----------------------------------------------------------------------*
rc = c_close(desc)
if (rc < 0) then
  AixCls = AixErr(ErrTxt)
  call SysWarnFileMsg('AixCls',FCtlBlk(nFile),'MSG: close',ErrTxt)
  call Abend()
end if
!----------------------------------------------------------------------*
! Update control block                                                 *
!----------------------------------------------------------------------*
CtlBlk(pHndle,nFile) = 0
CtlBlk(pDesc,nFile) = 0
CtlBlk(pWhere,nFile) = 0
CtlBlk(pStat,nFile) = 0
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixCls
