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
! rc=AixPWr(Handle,Buf,nBuf,iDisk)                                     *
!                                                                      *
! A buffer is written to a file associated with the file handle. The   *
! operation is synchronous, but can be done in parallel by multiple    *
! threads                                                              *
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
! Author:  Victor P. Vysotskiy                                         *
!          University of Lund, Sweden, 2012-2013                       *
! Written: 2012-2013                                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!                                                                      *
!***********************************************************************

function AixPWr(handle,Buf,nBuf,iDisk)

use Fast_IO, only: CtlBlk, eEof, eNtOpn, FCtlBlk, MxFile, pDesc, pHndle, ProfData, pWhere
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: AixPWr
integer(kind=iwp), intent(in) :: handle, Buf(*), nBuf, iDisk
integer(kind=iwp) :: desc, Lu, n, nFile, pDisk, rc
real(kind=wp) :: CPUA, CPUE, TIOA, TIOE
character(len=80) :: ErrTxt
character(len=*), parameter :: TheName = 'AixPWr'
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
! Entry to AixPWr                                                      *
!----------------------------------------------------------------------*
AixPWr = 0
rc = 0
!----------------------------------------------------------------------*
! Check if file is opened.                                             *
!----------------------------------------------------------------------*
n = 1
do
  if (CtlBlk(pHndle,n) == handle) exit
  n = n+1
  if (n > MxFile) then
    AixPWr = eNtOpn
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
  ProfData(7,Lu) = ProfData(7,Lu)+1
end if
!----------------------------------------------------------------------*
! Write to file                                                        *
!----------------------------------------------------------------------*
CtlBlk(pWhere,nFile) = pDisk+nBuf
if (nBuf > 0) rc = c_pwrite_wrapper(desc,Buf,nBuf,pDisk)
if (rc < 0) then
  call FASTIO('STATUS')
  AixPWr = AixErr(ErrTxt)
  call SysQuitFileMsg(_RC_IO_ERROR_WRITE_,TheName,FCtlBlk(nFile),'Premature abort while writing buffer to disk: ',ErrTxt)
else if (rc /= nBuf) then
  call FASTIO('STATUS')
  AixPWr = eEof
  call SysQuitFileMsg(_RC_IO_ERROR_WRITE_,TheName,FCtlBlk(nFile),'Premature abort while writing buffer to disk: ','Disk full? ')
end if
call Timing(CPUA,CPUE,TIOA,TIOE)
ProfData(1,Lu) = ProfData(1,Lu)+1
ProfData(2,Lu) = ProfData(2,Lu)+nBuf
ProfData(3,Lu) = ProfData(3,Lu)+TIOE
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

contains

function c_pwrite_wrapper(FileDescriptor,Buffer,nBytes,Offset)

  use, intrinsic :: iso_c_binding, only: c_loc

  integer(kind=iwp) :: c_pwrite_wrapper
  integer(kind=iwp), intent(in) :: FileDescriptor, nBytes, Offset
  integer(kind=iwp), intent(in), target :: Buffer(*)
  interface
    function c_pwrite(FileDescriptor,Buffer,nBytes,Offset) bind(C,name='c_pwrite_')
      use, intrinsic :: iso_c_binding, only: c_ptr
      use Definitions, only: MOLCAS_C_INT
      integer(kind=MOLCAS_C_INT) :: c_pwrite
      type(c_ptr), value :: Buffer
      integer(kind=MOLCAS_C_INT) :: FileDescriptor, nBytes, Offset
    end function c_pwrite
  end interface

  c_pwrite_wrapper = c_pwrite(FileDescriptor,c_loc(Buffer(1)),nBytes,Offset)

end function c_pwrite_wrapper

end function AixPWr
