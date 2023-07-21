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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Vec(iVec1,nVec,Buf,lBuf,nDim,iOpt)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: write (iOpt=1) or read (iOpt=2) "new" vectors to buffer.

use ChoMP2_dec, only: NowSym

implicit none
integer iVec1, nVec, lBuf, nDim, iOpt
real*8 Buf(lBuf)
#include "chomp2.fh"
character(len=3), parameter :: ThisNm = 'Vec'
character(len=10), parameter :: SecNam = 'ChoMP2_Vec'
integer :: iSym, iJob, lTot, iAdr
logical DoClose

iSym = NowSym
DoClose = .false.

if (iOpt == 1) then

  if (lUnit_F(iSym,2) < 1) then
    call ChoMP2_OpenF(1,2,iSym)
    DoClose = .true.
  end if

  iJob = 1
  lTot = nDim*nVec
  iAdr = nDim*(iVec1-1)+1
  call ddaFile(lUnit_F(iSym,2),iJob,Buf,lTot,iAdr)

else if (iOpt == 2) then

  if (lUnit_F(iSym,2) < 1) then
    call ChoMP2_OpenF(1,2,iSym)
    DoClose = .true.
  end if

  iJob = 2
  lTot = nDim*nVec
  iAdr = nDim*(iVec1-1)+1
  call ddaFile(lUnit_F(iSym,2),iJob,Buf,lTot,iAdr)

else

  write(6,*) SecNam,': illegal option: iOpt = ',iOpt
  call ChoMP2_Quit(SecNam,'illegal option',' ')

end if

if (DoClose) then
  call ChoMP2_OpenF(2,2,iSym)
  DoClose = .false.
end if

end subroutine ChoMP2_Vec
