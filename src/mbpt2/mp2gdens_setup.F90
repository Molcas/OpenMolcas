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

use MBPT2_Global, only: Density, DiaA, ip_Density, ip_DiaA, ip_Mp2Lagr, ip_WDensity, mAdDel, mAdFro, mAdOcc, mAdVir, Mp2Lagr, &
                        WDensity
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, iSym, l_Density, l_DiaA, l_Mp2Lagr, nExtT, nOccT
#include "corbinf.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
l_Density = 0
l_Mp2Lagr = 0
l_DiaA = 0
do iSym=1,nSym
  l_Density = l_Density+(nOrb(iSym)+nDel(iSym))*(nOrb(iSym)+nDel(iSym))
  l_Mp2Lagr = l_Mp2Lagr+(nFro(iSym)+nOcc(iSym))*(nExt(iSym)+nDel(iSym))
  l_DiaA = l_DiaA+(nFro(iSym)+nOcc(iSym))*(nExt(iSym)+nDel(iSym))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Density,l_Density,label='MP2Density')
call mma_allocate(WDensity,l_Density,label='MP2WDensity')
call mma_allocate(Mp2Lagr,l_Mp2Lagr,label='MP2Lagr')
call mma_allocate(DiaA,l_DiaA,label='MP2DiaA')

Density(:) = Zero
WDensity(:) = Zero
Mp2Lagr(:) = Zero
DiaA(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
ip_Density(1) = 1
ip_WDensity(1) = 1
ip_Mp2Lagr(1) = 1
ip_DiaA(1) = 1
do i=1,nSym-1
  ip_Density(i+1) = ip_Density(i)+(nOrb(i)+nDel(i))*(nOrb(i)+nDel(i))
  ip_WDensity(i+1) = ip_WDensity(i)+(nOrb(i)+nDel(i))*(nOrb(i)+nDel(i))
  ip_Mp2Lagr(i+1) = ip_Mp2Lagr(i)+(nFro(i)+nOcc(i))*(nExt(i)+nDel(i))
  ip_DiaA(i+1) = ip_DiaA(i)+(nFro(i)+nOcc(i))*(nExt(i)+nDel(i))
end do

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
