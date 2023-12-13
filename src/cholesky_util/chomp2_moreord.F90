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

subroutine ChoMP2_MOReOrd(CMO,COcc,CVir)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: reorder MOs,
!
!          CMO(alpha,i) -> COcc(i,alpha)
!          CMO(alpha,a) -> CVir(alpha,a)

use Cholesky, only: nBas, nSym
use ChoMP2, only: iAOVir, iT1AOT, nFro, nOcc, nVir
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(_OUT_) :: COcc(*), CVir(*)
integer(kind=iwp) :: i, iCount, iSym, jCount, kOff1, kOff2

iCount = 0
do iSym=1,nSym

  jCount = iCount+nBas(iSym)*nFro(iSym)

  do i=1,nOcc(iSym)
    kOff1 = jCount+nBas(iSym)*(i-1)+1
    kOff2 = iT1AOT(iSym,iSym)+i
    call dCopy_(nBas(iSym),CMO(kOff1),1,COcc(kOff2),nOcc(iSym))
  end do

  kOff1 = jCount+nBas(iSym)*nOcc(iSym)
  kOff2 = iAOVir(iSym,iSym)
  CVir(kOff2+1:kOff2+nBas(iSym)*nVir(iSym)) = CMO(kOff1+1:kOff1+nBas(iSym)*nVir(iSym))

  iCount = iCount+nBas(iSym)*nBas(iSym)

end do

end subroutine ChoMP2_MOReOrd
