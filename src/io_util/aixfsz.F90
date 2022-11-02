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
! Copyright (C) 2012,2013, Victor P. Vysotskiy                         *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                             A I X - I / O                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixFsz(Handle)                                                    *
!                                                                      *
! Return size of file.                                                 *
!                                                                      *
! Input:  Handle   - This is the unique file identifier associated     *
!                    with the file. It is created by AixOpn, and must  *
!                    be used on subsequent references to the file.     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Victor P. Vysotskiy                                         *
!          Lund University, Sweden                                     *
! Written: 2012-2013                                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!                                                                      *
!***********************************************************************

function AixFsz(handle)

use Fast_IO, only: CtlBlk, eNtOpn, FCtlBlk, pDesc, pHndle, MxFile
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: AixFsz
integer(kind=iwp), intent(in) :: handle
integer(kind=iwp) :: desc, n, nFile, rc
character(len=80) :: ErrTxt
integer(kind=iwp), external :: c_stat
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
end interface

!----------------------------------------------------------------------*
! Entry to AixFsz                                                      *
!----------------------------------------------------------------------*
AixFsz = 0
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
do
  if (CtlBlk(pHndle,n) == handle) exit
  n = n+1
  if (n > MxFile) then
    AixFsz = eNtOpn
    return
  end if
end do
nFile = n
desc = CtlBlk(pDesc,nFile)
!----------------------------------------------------------------------*
! Get file size                                                        *
!----------------------------------------------------------------------*
rc = c_stat(desc)
if (rc < 0) then
  AixFsz = AixErr(ErrTxt)
  call SysWarnFileMsg('AixFsz',FCtlBlk(nFile),'MSG: close',ErrTxt)
  call Abend()
end if
AixFsz = rc
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixFsz
