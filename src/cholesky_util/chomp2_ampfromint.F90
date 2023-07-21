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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_AmpFromInt(Col,nDim,iCol,nCol,EOcc,EVir)
!
! Thomas Bondo Pedersen, Dec. 2007.
!
! Purpose: scale integrals with orbital energies to get
!          (minus) MP2 amplitudes: (ai|bj)/[e(a)-e(i)+e(b)-e(j)].

use ChoMP2_dec, only: NowSym

implicit none
integer nDim, nCol
real*8 Col(nDim,nCol), EOcc(*), EVir(*)
integer iCol(nCol)
#include "chomp2.fh"
#include "cholesky.fh"
integer iSym, bj_, bj, b, iSymb, j, iSymj, iSymi, iSyma, i, ai0
integer a, ai
real*8 Ebj, DE
integer MulD2h, k, l
! Statement function
MulD2h(k,l) = ieor(k-1,l-1)+1

iSym = NowSym
do bj_=1,nCol
  bj = iCol(bj_)
  call ChoMP2_Col_Invai(bj,iSym,b,iSymb,j,iSymj)
  Ebj = EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
  do iSymi=1,nSym
    iSyma = MulD2h(iSymi,iSym)
    do i=1,nOcc(iSymi)
      ai0 = iT1am(iSyma,iSymi)+nVir(iSyma)*(i-1)
      do a=1,nVir(iSyma)
        ai = ai0+a
        DE = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+Ebj
        Col(ai,bj_) = Col(ai,bj_)/DE
      end do
    end do
  end do
end do

end subroutine ChoMP2_AmpFromInt
