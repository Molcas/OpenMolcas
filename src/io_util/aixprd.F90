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

integer function AixPRd(handle,Buf,nBuf,iDisk,iErrSkip)
implicit integer(a-z)
#include "SysDef.fh"
#include "blksize.fh"
#include "switch.fh"
#include "ctl.fh"
dimension Buf(*)
character*80 ErrTxt
character*16 TheName
#ifndef _OLD_IO_STAT_
#include "pfio.fh"
real*8 CPUA, CPUE, TIOA, TIOE
#endif
#include "warnings.fh"
data TheName/'AixPRd'/

!----------------------------------------------------------------------*
! Entry to AixPRd                                                      *
!----------------------------------------------------------------------*
AixPRd = 0
rc = 0
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
100 continue
if (CtlBlk(pHndle,n) /= handle) then
  n = n+1
  if (n > MxFile) then
    AixPRd = eNtOpn
    return
  end if
  Go To 100
end if
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
