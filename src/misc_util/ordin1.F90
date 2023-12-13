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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine OrdIn1(iOpt,Buf0,lBuf0,iBatch)
!***********************************************************************
!                                                                      *
!     Purpose: Read a buffer of ordered two electron integrals         *
!                                                                      *
!     Note:    This subroutine has internal buffers.                   *
!                                                                      *
!    Calling parameters:                                               *
!    Buf0   : contains on output the integrals                         *
!    lBuf0  : number of integrals to be transfered                     *
!    iOpt   : option code (iOpt=1:start reading at first integral)     *
!                         (iOpt=2:continue reading)                    *
!    rc     : return code                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.O. Widmark                                 *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
!----------------------------------------------------------------------*
!                                                                      *
!     This is a buffer exclusively used for I/O buffering              *
!     of the ordered 2-el. integrals file                              *
!                                                                      *
!     !!!     The current version uses one buffer only     !!!         *
!     !!! Double buffering and asynchronous I/O is preseen !!!         *
!                                                                      *
!----------------------------------------------------------------------*

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use TwoDat, only: AuxTwo, isDAdr, lStRec, lTop, TocTwo
use Definitions, only: wp, iwp, RtoB

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iOpt, iBatch
real(kind=wp), intent(_OUT_) :: Buf0(*)
integer(kind=iwp), intent(inout) :: lBuf0
integer(kind=iwp) :: iDisk1, isBuf0, isBuf1, jOpt, kOpt = 0, lBuf1, LuTwo, nByte, ncopy, nleft
character :: Buf1(8*lStRec) = ' '
real(kind=wp), external :: C2R8

!----------------------------------------------------------------------*
! Fetch the unit number, disk start address and pointers               *
!----------------------------------------------------------------------*
LuTwo = AuxTwo%Unt
iDisk1 = AuxTwo%DaDa
isBuf1 = AuxTwo%Upk8
lBuf1 = AuxTwo%lBf1
if (iOpt == 1) then
  !--------------------------------------------------------------------*
  ! If this is the first block of a symmetry batch                     *
  ! get the disk disk start address and load the buffer                *
  !--------------------------------------------------------------------*
  iDisk1 = TocTwo(isDAdr+iBatch-1)
  jOpt = 2
  call cdDAFILE(LuTwo,jOpt,Buf1,lStRec,iDisk1)
  !--------------------------------------------------------------------*
  ! Note:                                                              *
  ! If the records are organized in sequential order                   *
  ! (SORT3 in SEWARD is activated)                                     *
  ! deactivate the update of the disk address                          *
  !                                                                    *
  ! iDisk1=NINT(C2R8(Buf1(1)))                                         *
  !--------------------------------------------------------------------*
  lBuf1 = nint(C2R8(Buf1(17)))
  kOpt = nint(C2R8(Buf1(25)))
  isBuf1 = lTop*RtoB+1
end if
if (lBuf0 <= lBuf1) then
  !--------------------------------------------------------------------*
  ! If the number of requested integrals is smaller than               *
  ! the current buffer transfer data                                   *
  !--------------------------------------------------------------------*
  call cUPKR8(kOpt,lBuf0,nByte,Buf1(isBuf1),Buf0)
  isBuf1 = isBuf1+nByte
  lBuf1 = lBuf1-lBuf0
else
  !--------------------------------------------------------------------*
  ! If the number of requested integrals is larger than                *
  ! the current buffer first drain the current buffer and              *
  ! read as many subsequent buffers as needed                          *
  !--------------------------------------------------------------------*
  call cUPKR8(kOpt,lBuf1,nByte,Buf1(isBuf1),Buf0)
  isBuf0 = lBuf1+1
  nleft = lBuf0-lBuf1
  do while (nleft > 0)
    jOpt = 2
    call cdDAFILE(LuTwo,jOpt,Buf1,lStRec,iDisk1)
    !------------------------------------------------------------------*
    ! Note:                                                            *
    ! If the records are organized in sequential order                 *
    ! (SORT3 in SEWARD is activated)                                   *
    ! deactivate the update of the disk address                        *
    !                                                                  *
    ! iDisk1=NINT(C2R8(Buf1(1)))                                       *
    !------------------------------------------------------------------*
    lBuf1 = nint(C2R8(Buf1(17)))
    ncopy = min(nleft,lBuf1)
    kOpt = nint(C2R8(Buf1(25)))
    isBuf1 = lTop*RtoB+1
    call cUPKR8(kOpt,ncopy,nByte,Buf1(isBuf1),Buf0(isBuf0))
    isBuf0 = isBuf0+ncopy
    isBuf1 = lTop*RtoB+1+nByte
    lBuf1 = lBuf1-ncopy
    nleft = nleft-ncopy
  end do
end if
!----------------------------------------------------------------------*
! Update pointer to next disk address and integral to unpack           *
!----------------------------------------------------------------------*
AuxTwo%DaDa = iDisk1
AuxTwo%Upk8 = isBuf1
AuxTwo%lBf1 = lBuf1

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

!This is to allow type punning without an explicit interface
contains

subroutine cdDAFILE(Lu,iOpt,Buf,lBuf_,iDisk_)

  integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf_
  character, target, intent(inout) :: Buf(*)
  integer(kind=iwp), intent(inout) :: iDisk_
  real(kind=wp), pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf(1)),pBuf,[lBuf_])
  call dDAFILE(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine cdDAFILE

subroutine cUPKR8(iOpt,nData,nByte,InBuf,OutBuf)

  integer(kind=iwp), intent(in) :: iOpt
  integer(kind=iwp), intent(inout) :: nData
  integer(kind=iwp), intent(out) :: nByte
  character, target, intent(in) :: InBuf(*)
  real(kind=wp), intent(_OUT_) :: OutBuf(*)
  real(kind=wp), pointer :: dBuf(:)

  call c_f_pointer(c_loc(InBuf(1)),dBuf,[nData])
  call UPKR8(iOpt,nData,nByte,dBuf,OutBuf)
  nullify(dBuf)

end subroutine cUPKR8

end subroutine OrdIn1
