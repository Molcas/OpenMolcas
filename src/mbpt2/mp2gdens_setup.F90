!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine MP2gDens_setup()

use MBPT2_Global, only: Density, DiaA, mAdDel, mAdFro, mAdOcc, mAdVir, Mp2Lagr, WDensity
use Data_Structures, only: Allocate_DT
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iSym, nExtT, nOccT
#include "corbinf.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
call Allocate_DT(Density,nOrb+nDel,nOrb+nDel,nSym,label='MP2Density')
call Allocate_DT(WDensity,nOrb+nDel,nOrb+nDel,nSym,label='MP2WDensity')
call Allocate_DT(Mp2Lagr,nFro+nOcc,nExt+nDel,nSym,label='MP2Lagr')
call Allocate_DT(DiaA,nFro+nOcc,nExt+nDel,nSym,Label='MP2DiaA')

Density%A0(:) = Zero
WDensity%A0(:) = Zero
Mp2Lagr%A0(:) = Zero
DiaA%A0(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
mAdOcc(1) = 1
nOccT = nOcc(1)
do iSym=2,nSym
  mAdOcc(iSym) = mAdOcc(iSym-1)+nOcc(iSym-1)
  nOccT = nOccT+nOcc(iSym)
end do

mAdVir(1) = 1
nExtT = nExt(1)
do iSym=2,nSym
  mAdVir(iSym) = mAdVir(iSym-1)+nExt(iSym-1)
  nExtT = nExtT+nExt(iSym)
end do

mAdFro(1) = nOccT+1
do iSym=2,nSym
  mAdFro(iSym) = mAdFro(iSym-1)+nFro(iSym-1)
end do

mAdDel(1) = nExtT+1
do iSym=2,nSym
  mAdDel(iSym) = mAdDel(iSym-1)+nDel(iSym-1)
end do

return

end subroutine MP2gDens_setup
