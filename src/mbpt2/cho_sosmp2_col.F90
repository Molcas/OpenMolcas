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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_Col(Col,nDim,iCol,nCol,Buf,l_Buf)
! Francesco Aquilante, May 2007.
!
! Purpose: compute specified M(ai,bj)=(ai|bj)^2 columns.

use ChoMP2, only: OldVec

#include "implicit.fh"
real*8 Col(nDim,nCol), Buf(l_Buf)
integer iCol(nCol)
character*3 ThisNm
character*14 SecNam
parameter(SecNam='Cho_SOSmp2_Col',ThisNm='Col')
logical DoClose
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_dec.fh"
#include "WrkSpc.fh"

if ((nCol < 1) .or. (nDim < 1)) return

iSym = NowSym
if (nDim /= nT1am(iSym)) then
  write(6,*) SecNam,': inconsistent dimension. Expected: ',nT1am(iSym),'   Received: ',nDim
  write(6,*) SecNam,': symmetry from chomp2_dec.fh: ',iSym
  call ChoMP2_Quit(SecNam,'inconsistent dimension',' ')
end if

if (NumCho(iSym) < 1) then
  call Cho_dZero(Col,nDim*nCol)
  return
end if

irc = 0

if (InCore(iSym)) then  ! old vectors available in core

  Fac = 0.0d0
  call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,OldVec,NumCho(iSym),Buf,l_Buf,Fac,irc)
  if (irc /= 0) then
    write(6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
    call ChoMP2_Quit(SecNam,'ChoMP2_Col_Comp error','[1]')
  end if

else ! old vectors must be read on disk

  DoClose = .false.
  if (lUnit_F(iSym,1) < 1) then
    call ChoMP2_OpenF(1,1,iSym)
    DoClose = .true.
  end if

  call GetMem('MaxCol','Max ','Real',ipWrk,lWrk)

  if (l_Buf > lWrk) then ! use Buf as work space

    nVec = min(l_Buf/(nDim+1),NumCho(iSym))
    if (nVec < 1) then
      write(6,*) SecNam,': insufficient memory for batch!'
      call ChoMP2_Quit(SecNam,'insufficient memory','[1]')
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
      call ddaFile(lUnit_F(iSym,1),iOpt,Buf(1),lTot,iAdr)

      if (iBat == 1) then
        Fac = 0.0d0
      else
        Fac = 1.0d0
      end if

      lScr = l_Buf-lTot
      if (lWrk > lScr) then
        lWsav = lWrk
        call GetMem('ColScr','Allo','Real',ipWrk,lWrk)
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Buf(1),NumV,Work(ipWrk),lWrk,Fac,irc)
        call GetMem('ColScr','Free','Real',ipWrk,lWrk)
        lWrk = lWsav
      else
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Buf(1),NumV,Buf(1+lTot),lScr,Fac,irc)
      end if
      if (irc /= 0) then
        write(6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
        call ChoMP2_Quit(SecNam,'ChoMP2_Col_Comp error','[2]')
      end if

    end do

  else ! use Work as work space

    call GetMem('ColWrk','Allo','Real',ipWrk,lWrk)

    nVec = min(lWrk/nDim,NumCho(iSym))
    if (nVec < 1) then
      write(6,*) SecNam,': insufficient memory for batch!'
      call ChoMP2_Quit(SecNam,'insufficient memory','[2]')
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
      call ddaFile(lUnit_F(iSym,1),iOpt,Work(ipWrk),lTot,iAdr)

      if (iBat == 1) then
        Fac = 0.0d0
      else
        Fac = 1.0d0
      end if

      lScr = lWrk-lTot
      if (l_Buf > lScr) then
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Work(ipWrk),NumV,Buf(1),l_Buf,Fac,irc)
      else
        call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Work(ipWrk),NumV,Work(ipWrk+lTot),lScr,Fac,irc)
      end if
      if (irc /= 0) then
        write(6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
        call ChoMP2_Quit(SecNam,'ChoMP2_Col_Comp error','[3]')
      end if

    end do

    call GetMem('ColWrk','Free','Real',ipWrk,lWrk)

  end if

  if (DoClose) then
    call ChoMP2_OpenF(2,1,iSym)
    DoClose = .false.
  end if

end if

! Squaring each element of the integral columns
! ---------------------------------------------
do jCol=1,nCol
  do ia=1,nDim
    Col(ia,jCol) = Col(ia,jCol)**2
  end do
end do

end subroutine Cho_SOSmp2_Col
