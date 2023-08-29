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

subroutine ChoMP2_Col(Col,nDim,iCol,nCol,Buf,l_Buf)
!
! Thomas Bondo Pedersen, Dec. 2007.
!
! Purpose: compute specified (ai|bj) or MP2 amplitude columns.

use ChoMP2, only: EOcc, EVir, iOption_MP2CD, NowSym, nT1am
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDim, nCol, iCol(nCol), l_Buf
real(kind=wp), intent(out) :: Col(nDim,nCol), Buf(l_Buf)
integer(kind=iwp) :: iSym
character(len=*), parameter :: SecNam = 'ChoMP2_Col'

! Check input.
! ------------

if ((nCol < 1) .or. (nDim < 1)) return
iSym = NowSym
if (nDim /= nT1am(iSym)) then
  write(u6,*) SecNam,': inconsistent dimension. Expected: ',nT1am(iSym),'   Received: ',nDim
  write(u6,*) SecNam,': symmetry from Module chomp2_dec: ',iSym
  call SysAbendMsg(SecNam,'inconsistent dimension',' ')
end if

! Calculate (ai|bj) integrals.
! ----------------------------

call ChoMP2_IntCol(Col,nDim,iCol,nCol,Buf,l_Buf)

! Postprocess integrals.
! ----------------------

if (iOption_MP2CD == 2) call ChoMP2_AmpFromInt(Col,nDim,iCol,nCol,EOcc,EVir) ! generate amplitudes

end subroutine ChoMP2_Col
