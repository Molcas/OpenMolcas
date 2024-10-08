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
! Copyright (C) 1998, Roland Lindh                                     *
!***********************************************************************

subroutine IniBuf(nDisc,nCore)
!***********************************************************************
!                                                                      *
!  Object: Initiate I/O buffer for semi-direct SCF                     *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemical Physics,                 *
!             University of Lund, Sweden. October '98                  *
!***********************************************************************

use IOBUF, only: Buffer, DiskMx_Byte, InCore, lBuf, LuTmp, nBuf, OnDisk
use stdalloc, only: mma_allocate, mma_maxDBLE
use Constants, only: Two, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDisc
integer(kind=iwp), intent(inout) :: nCore
integer(kind=iwp) :: MaxMem, MemMin_Seward, MemReq
integer(kind=iwp), external :: AllocDisk

! Open file for semi-direct implementation
! nDisc in units of MByte
! nCore in units of kByte

! The maximum number of bytes on disk. The file size limit times
! the number of multi files.

DiskMx_Byte = real(AllocDisk(),kind=wp)*Ten*Two**20

nBuf = -99
if ((nDisc == 0) .and. (nCore == 0)) then
  OnDisk = .false.
  Incore = .false.
else if (nDisc*1024 > nCore) then
  OnDisk = .true.
  Incore = .false.
  LuTMp = 32
  call EAFOpen(LuTmp,'SMDINT  ')
  nBuf = 2
else
  OnDisk = .false.
  InCore = .true.
  nBuf = 1
end if

! Adjust buffer size and allocate memory for the buffer.

if (OnDisk .or. InCore) then
  MemMin_Seward = 1024**2 ! Real
  call mma_maxDBLE(MaxMem)
  ! lBuf in units of Real per buffer
  lBuf = (1024*nCore)/(8*nBuf)
  if (InCore) then
    MemReq = MemMin_Seward+lBuf*nBuf
    !write(u6,*) 'MemReq,MaxMem=',MemReq,MaxMem,'  lbuf=',lbuf
    if (MemReq > MaxMem) then
      lBuf = (MaxMem-MemMin_Seward)/nBuf
      if (lBuf < 0) then
        nCore = (((MaxMem*3)/4)*8)/1024
      else
        nCore = (lBuf*8)/1024
      end if
    else
      nCore = (lBuf*8)/1024
    end if
    nCore = ((nCore+7)/8)*8
    lBuf = (1024*nCore)/(8*nBuf)
  end if
  !write(u6,*) 'OnDisk=',OnDisk
  !write(u6,*) 'Incore=',Incore
  !write(u6,*) 'nBuf=',nBuf
  !write(u6,*) 'IniBuf: nDisc=',nDisc,'MByte'
  !write(u6,*) 'nCore=',nCore,'kByte'
  !write(u6,*) 'lBuf=',lBuf

  ! Allocate I/O Buffer
  call mma_allocate(Buffer,lBuf,nBuf,Label='Buffer')
end if

end subroutine IniBuf
