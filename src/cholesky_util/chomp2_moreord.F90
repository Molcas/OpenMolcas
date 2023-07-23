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

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: COcc(*), CVir(*), CMO(*)
#include "cholesky.fh"
#include "chomp2.fh"
#include "choorb.fh"
integer(kind=iwp) :: i, iCount, iSym, jCount, kOff1, kOff2

iCount = 0
do iSym=1,nSym

  jCount = iCount+nBas(iSym)*nFro(iSym)

  do i=1,nOcc(iSym)
    kOff1 = jCount+nBas(iSym)*(i-1)+1
    kOff2 = iT1AOT(iSym,iSym)+i
    call dCopy_(nBas(iSym),CMO(kOff1),1,COcc(kOff2),nOcc(iSym))
  end do

  kOff1 = jCount+nBas(iSym)*nOcc(iSym)+1
  kOff2 = iAOVir(iSym,iSym)+1
  call dCopy_(nBas(isym)*nVir(iSym),CMO(kOff1),1,CVir(kOff2),1)

  iCount = iCount+nBas(iSym)*nBas(iSym)

end do

end subroutine ChoMP2_MOReOrd
