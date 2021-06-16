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
!               2012,2013, Victor P. Vysotskiy                         *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                             A I X - I / O                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixRd(Handle,Buf,nBuf,iDisk,iErrSkip)                             *
!                                                                      *
! A buffer is read from a file associated with the file handle. The    *
! operation is asynchronous, and must be followed by a call to AixWt   *
! to ensure that data is in memory.                                    *
!                                                                      *
! Input:  Handle   - This is the unique file identifier associated     *
!                    with the file. It is created by AixOpn, and must  *
!                    be used on subsequent references to the file.     *
!         Buf      - The buffer that is to be written to disk.         *
!         nBuf     - Length of the buffer in words.                    *
!         iDisk    - External disk address.                            *
!         iErrSkip - if 0 : stop in the error case                     *
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
!          V.P. Vysotskiy,University of Lund, Sweden, 2012-2013        *
!                                                                      *
!***********************************************************************

function AixRd(handle,Buf,nBuf,iDisk,iErrSkip)

#include "intent.fh"

use Fast_IO, only: CtlBlk, eEof, eInErr, eNtOpn, FCtlBlk, MxFile, pDesc, pHndle, ProfData, pWhere
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: AixRd
integer(kind=iwp), intent(in) :: handle, nBuf, iErrSkip
integer(kind=iwp), intent(_OUT_) :: Buf(*)
integer(kind=iwp), intent(inout) :: iDisk
integer(kind=iwp) :: desc, Lu, n, nFile, pDisk, rc
real(kind=wp) :: CPUA, CPUE, TIOA, TIOE
character(len=80) :: ErrTxt
character(len=64) :: Temp
character(len=*), parameter :: TheName = 'AixRd'
integer(kind=iwp), external :: AixErr, c_lseek, c_read
#include "warnings.h"

!----------------------------------------------------------------------*
! Entry to AixRd                                                       *
!----------------------------------------------------------------------*
AixRd = 0
Temp = 'Premature abort while reading buffer from disk'
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
do
  if (CtlBlk(pHndle,n) == handle) exit
  n = n+1
  if (n > MxFile) then
    AixRd = eNtOpn
    return
  end if
end do
nFile = n
desc = CtlBlk(pDesc,nFile)
call FSCB2UNIT(handle,Lu)
call Timing(CPUA,CPUE,TIOA,TIOE)
!----------------------------------------------------------------------*
! Position file pointer                                                *
!----------------------------------------------------------------------*
pDisk = iDisk
if (CtlBlk(pWhere,nFile) /= pDisk) then
  rc = c_lseek(desc,pDisk)
  ProfData(8,Lu) = ProfData(8,Lu)+1
  if (rc < 0) then
    if (iErrSkip == 1) then
      AixRd = 99
      return
    end if
    call FASTIO('STATUS')
    AixRd = AixErr(ErrTxt)
    call SysWarnFileMsg(TheName,FCtlBlk(nFile),'MSG: seek',ErrTxt)
    call SysCondMsg('rc < 0',rc,'<',0)
  else if (rc /= pDisk) then
    if (iErrSkip == 1) then
      AixRd = 99
      return
    end if
    call FASTIO('STATUS')
    AixRd = eInErr
    call SysWarnFileMsg(TheName,FCtlBlk(nFile),'MSG: seek',' ')
    call SysCondMsg('rc != pDisk',rc,'!=',pDisk)
  end if
end if
CtlBlk(pWhere,nFile) = pDisk
!----------------------------------------------------------------------*
! Read from file                                                       *
!----------------------------------------------------------------------*
rc = c_read(desc,Buf,nBuf)
if (rc < 0) then
  if (iErrSkip == 1) then
    AixRd = 99
    return
  end if
  call FASTIO('STATUS')
  AixRd = AixErr(ErrTxt)

  call SysQuitFileMsg(_RC_IO_ERROR_READ_,TheName,FCtlBlk(nFile),Temp,ErrTxt)

else if (rc /= nBuf) then
  if (iErrSkip == 1) then
    AixRd = 99
    return
  end if
  call FASTIO('STATUS')
  AixRd = eEof
  call SysQuitFileMsg(_RC_IO_ERROR_READ_,TheName,FCtlBlk(nFile),Temp,'\nEnd of file reached ')
end if
CtlBlk(pWhere,nFile) = CtlBlk(pWhere,nFile)+nBuf
iDisk = iDisk+nBuf
call Timing(CPUA,CPUE,TIOA,TIOE)
ProfData(4,Lu) = ProfData(4,Lu)+1
ProfData(5,Lu) = ProfData(5,Lu)+nBuf
ProfData(6,Lu) = ProfData(6,Lu)+TIOE
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixRd
