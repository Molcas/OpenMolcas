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

integer function AixCls(handle)

implicit integer(a-z)
#include "switch.fh"
#include "ctl.fh"
character*80 ErrTxt

!----------------------------------------------------------------------*
! Entry to AixCls                                                      *
!----------------------------------------------------------------------*
AixCls = 0
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
100 continue
if (CtlBlk(pHndle,n) /= handle) then
  n = n+1
  if (n > MxFile) then
    AixCls = eNtOpn
    return
  end if
  Go To 100
end if
nFile = n
desc = CtlBlk(pDesc,nFile)
!----------------------------------------------------------------------*
! Close file                                                           *
!----------------------------------------------------------------------*
rc = c_close(desc)
if (rc < 0) then
  AixCls = AixErr(ErrTxt)
  call SysAbendFileMsg('AixCls',FCtlBlk(nFile),'MSG: close',ErrTxt)
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
