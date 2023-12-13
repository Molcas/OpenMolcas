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
! Copyright (C) 2004,2007, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine ChoMP2_IntCol(Col,nDim,iCol,nCol,Buf,l_Buf)
!
! Thomas Bondo Pedersen, Dec. 2004.
! Renamed (from ChoMP2_Col), Thomas Bondo Pedersen, Dec. 2007.
!
! Purpose: compute specified (ai|bj) columns.

use Cholesky, only: NumCho
use ChoMP2, only: InCore, lUnit_F, NowSym, OldVec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDim, nCol, iCol(nCol), l_Buf
real(kind=wp), intent(out) :: Col(nDim,nCol), Buf(l_Buf)
integer(kind=iwp) :: iAdr, iBat, iOpt, irc, iSym, iVec1, lScr, lTot, lWrk, lWsav, nBat, NumV, nVec
real(kind=wp) :: Fac
logical(kind=iwp) :: DoClose
real(kind=wp), allocatable :: Wrk(:)
character(len=*), parameter :: SecNam = 'ChoMP2_IntCol'

iSym = NowSym
if (NumCho(iSym) < 1) then
  Col(:,:) = Zero
  return
end if

irc = 0

if (InCore(iSym)) then  ! old vectors available in core

  Fac = Zero
  call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,OldVec,NumCho(iSym),Buf,l_Buf,Fac,irc)
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
    call SysAbendMsg(SecNam,'ChoMP2_Col_Comp error','[1]')
  end if

else ! old vectors must be read on disk

  DoClose = .false.
  if (lUnit_F(iSym,1) < 1) then
    call ChoMP2_OpenF(1,1,iSym)
    DoClose = .true.
  end if

  call mma_maxDBLE(lWrk)

  if (l_Buf > lWrk) then ! use Buf as work space

    nVec = min(l_Buf/(nDim+1),NumCho(iSym))
    if (nVec < 1) then
      write(u6,*) SecNam,': insufficient memory for batch!'
      call SysAbendMsg(SecNam,'insufficient memory','[1]')
      nBat = 0
    else
      nBat = (NumCho(iSym)-1)/nVec+1
    end if

    do iBat=1,nBat

      if (iBat == nBat) then
        NumV = NumCho(iSym)-nVec*(nBat-1)
      else
        NumV = nVec
      end if
      iVec1 = nVec*(iBat-1)+1

      iOpt = 2
      lTot = nDim*NumV
      iAdr = nDim*(iVec1-1)+1
      call ddaFile(lUnit_F(iSym,1),iOpt,Buf,lTot,iAdr)

      if (iBat == 1) then
        Fac = Zero
      else
        Fac = One
      end if

      lScr = l_Buf-lTot
      if (lWrk > lScr) then
        lWsav = lWrk
        call mma_allocate(Wrk,lWrk,Label='Wrk')
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Buf(1),NumV,Wrk,lWrk,Fac,irc)
        call mma_deallocate(Wrk)
        lWrk = lWsav
      else
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Buf(1),NumV,Buf(1+lTot),lScr,Fac,irc)
      end if
      if (irc /= 0) then
        write(u6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
        call SysAbendMsg(SecNam,'ChoMP2_Col_Comp error','[2]')
      end if

    end do

  else ! use Work as work space

    call mma_allocate(Wrk,lWrk,Label='Wrk')

    nVec = min(lWrk/nDim,NumCho(iSym))
    if (nVec < 1) then
      write(u6,*) SecNam,': insufficient memory for batch!'
      call SysAbendMsg(SecNam,'insufficient memory','[2]')
      nBat = 0
    else
      nBat = (NumCho(iSym)-1)/nVec+1
    end if

    do iBat=1,nBat

      if (iBat == nBat) then
        NumV = NumCho(iSym)-nVec*(nBat-1)
      else
        NumV = nVec
      end if
      iVec1 = nVec*(iBat-1)+1

      iOpt = 2
      lTot = nDim*NumV
      iAdr = nDim*(iVec1-1)+1
      call ddaFile(lUnit_F(iSym,1),iOpt,Wrk,lTot,iAdr)

      if (iBat == 1) then
        Fac = Zero
      else
        Fac = ONe
      end if

      lScr = lWrk-lTot
      if (l_Buf > lScr) then
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Wrk,NumV,Buf(1),l_Buf,Fac,irc)
      else
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Wrk,NumV,Wrk(1+lTot),lScr,Fac,irc)
      end if
      if (irc /= 0) then
        write(u6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
        call SysAbendMsg(SecNam,'ChoMP2_Col_Comp error','[3]')
      end if

    end do

    call mma_deallocate(Wrk)

  end if

  if (DoClose) then
    call ChoMP2_OpenF(2,1,iSym)
    DoClose = .false.
  end if

end if

end subroutine ChoMP2_IntCol
