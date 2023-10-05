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

use Symmetry_Info, only: Mul
use Cholesky, only: nSym
use ChoMP2, only: iOcc, iT1am, iVir, nOcc, NowSym, nVir
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nCol, iCol(nCol)
real(kind=wp), intent(inout) :: Col(nDim,nCol)
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
integer(kind=iwp) :: a, ai, ai0, b, bj, bj_, i, iSym, iSyma, iSymb, iSymi, iSymj, j
real(kind=wp) :: DE, Ebj

iSym = NowSym
do bj_=1,nCol
  bj = iCol(bj_)
  call ChoMP2_Col_Invai(bj,iSym,b,iSymb,j,iSymj)
  Ebj = EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
  do iSymi=1,nSym
    iSyma = Mul(iSymi,iSym)
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
