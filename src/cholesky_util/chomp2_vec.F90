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

use ChoMP2, only: lUnit_F, NowSym
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iVec1, nVec, lBuf, nDim, iOpt
real(kind=wp), intent(inout) :: Buf(lBuf)
integer(kind=iwp) :: iAdr, iJob, iSym, lTot
logical(kind=iwp) :: DoClose
character(len=*), parameter :: SecNam = 'ChoMP2_Vec'

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

  write(u6,*) SecNam,': illegal option: iOpt = ',iOpt
  call SysAbendMsg(SecNam,'illegal option',' ')

end if

if (DoClose) then
  call ChoMP2_OpenF(2,2,iSym)
  DoClose = .false.
end if

end subroutine ChoMP2_Vec
