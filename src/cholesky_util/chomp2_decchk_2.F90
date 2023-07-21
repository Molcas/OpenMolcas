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
! Copyright (C) 2008, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_DecChk_2(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
!
! Thomas Bondo Pedersen, Jan. 2008.
!
! Purpose: check MP2 decomposition of the MP2 amplitudes
!          (sym. block iSym). The columns of the matrix are
!          compared nCol columns at a time. This
!          implies that the memory requirement of this routine
!          should be limited to approximately that of the
!          decomposition itself. Note, however, that since all
!          integrals are computed, this routine will consume
!          significantly more CPU time.
!          Files are assumed open.
!          On exit,
!          ErrStat(1) = min error
!          ErrStat(2) = max error
!          ErrStat(3) = rms error

use ChoMP2, only: OldVec
use ChoMP2_dec, only: EOcc, EVir, Incore

implicit real*8(a-h,o-z)
real*8 Col(nDim,nCol), Wrk(lWrk), ErrStat(3)
#include "cholesky.fh"
#include "chomp2.fh"
external ddot_
integer a, b
character*8 ThisNm
character*15 SecNam
parameter(SecNam='ChoMP2_DecChk_2',ThisNm='DecChk_2')
! Statement function
MulD2h(k,l) = ieor(k-1,l-1)+1

irc = 0

! Check dimensions.
! -----------------

if ((nDim < 1) .or. (nCol < 1)) return
if (nDim /= nT1am(iSym)) then
  irc = -1
  Go To 1 ! exit
end if

! Initialize.
! -----------

ErrStat(1) = 9.9d15
ErrStat(2) = -9.9d15
ErrStat(3) = 0.0d0

! Set up batching over columns of the (ai|bj) matrix.
! ---------------------------------------------------

NumCol = min(nCol,nT1am(iSym))
nBatCol = (nT1am(iSym)-1)/NumCol+1

! Start column batch loop.
! ------------------------

Nai = nDim
do iBatCol=1,nBatCol

  ! Set batch info.
  ! ---------------

  if (iBatCol == nBatCol) then
    Nbj = nT1am(iSym)-NumCol*(nBatCol-1)
  else
    Nbj = NumCol
  end if
  ibj1 = NumCol*(iBatCol-1)+1

  ! Compute amplitudes from "old" vectors.
  ! --------------------------------------

  if (InCore(iSym)) then
    call DGEMM_('N','T',Nai,Nbj,NumCho(iSym),1.0d0,OldVec,Nai,OldVec(ibj1),Nai,0.0d0,Col,Nai)
  else
    lU = lUnit_F(iSym,1)
    NumVec = NumCho(iSym)
    Fac = 0.0d0
    call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,Wrk,lWrk,Fac)
    if (irc /= 0) then
      write(6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,' [2]'
      irc = 2
      Go To 1
    end if
  end if

  ibj0 = ibj1-1
  do iCol=1,Nbj
    ibj = ibj0+iCol
    call ChoMP2_Col_Invai(ibj,iSym,b,iSymb,j,iSymj)
    Ebj = EVir(iVir(iSymb)+b)-Eocc(iOcc(iSymj)+j)
    do iSymi=1,nSym
      iSyma = MulD2h(iSymi,iSym)
      do i=1,nOcc(iSymi)
        iai0 = iT1am(iSyma,iSymi)+nVir(iSyma)*(i-1)
        do a=1,nVir(iSyma)
          iai = iai0+a
          DE = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+Ebj
          Col(iai,iCol) = Col(iai,iCol)/DE
        end do
      end do
    end do
  end do

  ! Compute amplitudes from "new" vectors and subtract "old".
  ! ---------------------------------------------------------

  lU = lUnit_F(iSym,2)
  NumVec = nMP2Vec(iSym)
  Fac = -1.0d0
  call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,Wrk,lWrk,Fac)
  if (irc /= 0) then
    write(6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,' [1]'
    irc = 1
    Go To 1
  end if

  ! Compute error stats.
  ! --------------------

  do kbj=1,Nbj
    do kai=1,Nai
      ErrStat(1) = min(ErrStat(1),Col(kai,kbj))
      ErrStat(2) = max(ErrStat(2),Col(kai,kbj))
    end do
  end do
  ErrStat(3) = ErrStat(3)+dDot_(Nai*Nbj,Col,1,Col,1)

end do

! Compute rms error.
! ------------------

xdim = dble(Nai)*dble(Nai)
ErrStat(3) = sqrt(ErrStat(3)/xdim)

1 continue

end subroutine ChoMP2_DecChk_2
