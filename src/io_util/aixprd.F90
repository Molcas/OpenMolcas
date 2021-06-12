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
!                     Thread safe                                      *
!                     POSIX - I/O                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixPRd(Handle,Buf,nBuf,iDisk,iErrSkip)                            *
!                                                                      *
! A buffer is read from a file associated with the file handle. The    *
! operation is synchronous, but can be done in parallel by multiple    *
! threads                                                              *
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
! Author:  Victor P. Vysotskiy                                         *
!          University of Lund, Sweden, 2012-2013                       *
! Written: 2012-2013                                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!                                                                      *
!***********************************************************************

function AixPRd(handle,Buf,nBuf,iDisk,iErrSkip)

#include "intent.fh"

use Definitions, only: iwp
#ifndef _OLD_IO_STAT_
use Definitions, only: wp
#endif

implicit none
integer(kind=iwp) :: AixPRd
integer(kind=iwp), intent(in) :: handle, nBuf, iDisk, iErrSkip
integer(kind=iwp), intent(_OUT_) :: Buf(*)
integer(kind=iwp) :: desc, Lu, n, nFile, pDisk, rc
character(len=80) :: ErrTxt
character(len=6), parameter :: TheName = 'AixPRd'
integer(kind=iwp), external :: AixErr, c_pread
#include "switch.fh"
#include "ctl.fh"
#include "warnings.fh"
#ifndef _OLD_IO_STAT_
real(kind=wp) :: CPUA, CPUE, TIOA, TIOE
#include "pfio.fh"
#endif

!----------------------------------------------------------------------*
! Entry to AixPRd                                                      *
!----------------------------------------------------------------------*
AixPRd = 0
rc = 0
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
do
  if (CtlBlk(pHndle,n) == handle) exit
  n = n+1
  if (n > MxFile) then
    AixPRd = eNtOpn
    return
  end if
end do
nFile = n
desc = CtlBlk(pDesc,nFile)
#ifndef _OLD_IO_STAT_
call FSCB2UNIT(handle,Lu)
call Timing(CPUA,CPUE,TIOA,TIOE)
#endif
!----------------------------------------------------------------------*
! Position file pointer                                                *
!----------------------------------------------------------------------*
pDisk = pHeadOffset+iDisk
#ifndef _OLD_IO_STAT_
if (CtlBlk(pWhere,nFile) /= pDisk) then
  ProfData(8,Lu) = ProfData(8,Lu)+1
end if
#endif
!----------------------------------------------------------------------*
! Read from file                                                       *
!----------------------------------------------------------------------*
CtlBlk(pWhere,nFile) = pDisk+nBuf
if (nBuf > 0) rc = c_pread(desc,Buf,nBuf,pDisk)
if (rc < 0) then
  if (iErrSkip == 1) then
    AixPRd = 99
    return
  end if
  call FASTIO('STATUS')
  AixPRd = AixErr(ErrTxt)
  call SysQuitFileMsg(_RC_IO_ERROR_READ_,TheName,FCtlBlk(nFile),'Premature abort while reading buffer from disk',ErrTxt)

else if (rc /= nBuf) then
  if (iErrSkip == 1) then
    AixPRd = 99
    return
  end if
  call FASTIO('STATUS')
  AixPRd = eEof
  call SysQuitFileMsg(_RC_IO_ERROR_READ_,TheName,FCtlBlk(nFile),'Premature abort while reading buffer from disk:','\n End of file reached ')
end if
#ifndef _OLD_IO_STAT_
call Timing(CPUA,CPUE,TIOA,TIOE)
ProfData(4,Lu) = ProfData(4,Lu)+1
ProfData(5,Lu) = ProfData(5,Lu)+nBuf
ProfData(6,Lu) = ProfData(6,Lu)+TIOE
#endif
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixPRd
