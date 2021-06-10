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
!               2012, Victor P. Vysotskiy                              *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                             A I X - I / O                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixWr(Handle,Buf,nBuf,iDisk)                                      *
!                                                                      *
! A buffer is written to a file associated with the file handle. The   *
! operation is asynchronous, and must be followed by a call to AixWt   *
! to assure data integrity.                                            *
!                                                                      *
! Input:  Handle   - This is the unique file identifier associated     *
!                    with the file. It is created by AixOpn, and must  *
!                    be used on subsequent references to the file.     *
!         Buf      - The buffer that is to be written to disk.         *
!         nBuf     - Length of the buffer in words.                    *
!         iDisk    - External disk address.                            *
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
!     V.P. Vysotskiy,University of Lund, Sweden, 2012                  *
!                                                                      *
!***********************************************************************

integer function AixWr(handle,Buf,nBuf,iDisk)
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
data TheName/'AixWr'/

!----------------------------------------------------------------------*
! Entry to AixWr                                                       *
!----------------------------------------------------------------------*
AixWr = 0
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
100 continue
if (CtlBlk(pHndle,n) /= handle) then
  n = n+1
  if (n > MxFile) then
    AixWr = eNtOpn
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
if (CtlBlk(pWhere,nFile) /= pDisk) then
  rc = c_lseek(desc,pDisk)
# ifndef _OLD_IO_STAT_
  ProfData(7,Lu) = ProfData(7,Lu)+1
# endif
  if (rc < 0) then
    call FASTIO('STATUS')
    AixWr = AixErr(ErrTxt)
    call SysWarnFileMsg(TheName,FCtlBlk(nFile),'MSG: seek',ErrTxt)
    call SysCondMsg('rc < 0',rc,'<',0)
  else if (rc /= pDisk) then
    call FASTIO('STATUS')
    AixWr = eInErr
    call SysWarnFileMsg(TheName,FCtlBlk(nFile),'MSG: seek',' ')
    call SysCondMsg('rc != pDisk',rc,'!=',pDisk)
  end if
end if
CtlBlk(pWhere,nFile) = pDisk
!----------------------------------------------------------------------*
! Write to file                                                        *
!----------------------------------------------------------------------*
rc = c_write(desc,Buf,nBuf)
if (rc < 0) then
  call FASTIO('STATUS')
  AixWr = AixErr(ErrTxt)
  call SysQuitFileMsg(_RC_IO_ERROR_WRITE_,TheName,FCtlBlk(nFile),'Premature abort while writing buffer to disk:',ErrTxt)
else if (rc /= nBuf) then
  call FASTIO('STATUS')
  AixWr = eEof
  call SysQuitFileMsg(_RC_IO_ERROR_WRITE_,TheName,FCtlBlk(nFile),'Premature abort while writing buffer to disk:','Disk full? ')
end if
CtlBlk(pWhere,nFile) = CtlBlk(pWhere,nFile)+nBuf
iDisk = iDisk+nBuf
#ifndef _OLD_IO_STAT_
call Timing(CPUA,CPUE,TIOA,TIOE)
ProfData(1,Lu) = ProfData(1,Lu)+1
ProfData(2,Lu) = ProfData(2,Lu)+nBuf
ProfData(3,Lu) = ProfData(3,Lu)+TIOE
#endif
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixWr
