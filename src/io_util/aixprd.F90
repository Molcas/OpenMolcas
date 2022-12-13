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

use Fast_IO, only: CtlBlk, eEof, eNtOpn, FCtlBlk, MxFile, pDesc, pHndle, ProfData, pWhere
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: AixPRd
integer(kind=iwp), intent(in) :: handle, nBuf, iDisk, iErrSkip
integer(kind=iwp), intent(_OUT_) :: Buf(*)
integer(kind=iwp) :: desc, Lu, n, nFile, pDisk, rc
real(kind=wp) :: CPUA, CPUE, TIOA, TIOE
character(len=80) :: ErrTxt
character(len=*), parameter :: TheName = 'AixPRd'
#include "warnings.h"
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
end interface

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
call FSCB2UNIT(handle,Lu)
call Timing(CPUA,CPUE,TIOA,TIOE)
!----------------------------------------------------------------------*
! Position file pointer                                                *
!----------------------------------------------------------------------*
pDisk = iDisk
if (CtlBlk(pWhere,nFile) /= pDisk) then
  ProfData(8,Lu) = ProfData(8,Lu)+1
end if
!----------------------------------------------------------------------*
! Read from file                                                       *
!----------------------------------------------------------------------*
CtlBlk(pWhere,nFile) = pDisk+nBuf
if (nBuf > 0) rc = c_pread_wrapper(desc,Buf,nBuf,pDisk)
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
  call SysQuitFileMsg(_RC_IO_ERROR_READ_,TheName,FCtlBlk(nFile),'Premature abort while reading buffer from disk:', &
                      '\n End of file reached ')
end if
call Timing(CPUA,CPUE,TIOA,TIOE)
ProfData(4,Lu) = ProfData(4,Lu)+1
ProfData(5,Lu) = ProfData(5,Lu)+nBuf
ProfData(6,Lu) = ProfData(6,Lu)+TIOE
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

contains

function c_pread_wrapper(FileDescriptor,Buffer,nBytes,Offset)

  use, intrinsic :: iso_c_binding, only: c_loc

  integer(kind=iwp) :: c_pread_wrapper
  integer(kind=iwp), intent(in) :: FileDescriptor, nBytes, Offset
  integer(kind=iwp), intent(_OUT_), target :: Buffer(*)
  interface
    function c_pread(FileDescriptor,Buffer,nBytes,Offset) bind(C,name='c_pread_')
      use, intrinsic :: iso_c_binding, only: c_ptr
      use Definitions, only: MOLCAS_C_INT
      integer(kind=MOLCAS_C_INT) :: c_pread
      type(c_ptr), value :: Buffer
      integer(kind=MOLCAS_C_INT) :: FileDescriptor, nBytes, Offset
    end function c_pread
  end interface

  c_pread_wrapper = c_pread(FileDescriptor,c_loc(Buffer(1)),nBytes,Offset)

end function c_pread_wrapper

end function AixPRd
