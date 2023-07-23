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
! Copyright (C) 2008, Jonas Bostrom                                    *
!***********************************************************************

subroutine Construct_WDensII(EOcc,EVir,EFro,EDel)
!
! Jonas Bostrom, October 2008
!
! Purpose: Construct the piece of the energy-weighted density
!          usually labeled II.

use Constants, only: Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: EOcc(*), EVir(*), EFro(*), EDel(*)
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iA, iB, iI, iJ, iSym
real(kind=wp) :: E_a, E_b, E_i, E_j
! Statement functions
integer(kind=iwp) :: iDensActOcc, iWDensActOcc, iDensVactVall, iWDensVactVall, iDensVallOcc, iWDensVallOcc, i, j, k
iDensActOcc(i,j,k) = ip_Density(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)-1)
iWDensActOcc(i,j,k) = ip_WDensity(k)+j-1+(nOrb(k)+nDel(k))*(i-1+nFro(k))
iDensVactVall(i,j,k) = ip_Density(k)+j-1+nFro(k)+nOcc(k)+(nOrb(k)+nDel(k))*(i-1+nFro(k)+nOcc(k))
iWDensVactVall(i,j,k) = ip_WDensity(k)+j-1+nFro(k)+nOcc(k)+(nOrb(k)+nDel(k))*(i-1+nFro(k)+nOcc(k))
iDensVallOcc(i,j,k) = ip_Density(k)+j-1+(nOrb(k)+nDel(k))*(i-1+nFro(k)+nOcc(k))
iWDensVallOcc(i,j,k) = ip_WDensity(k)+j-1+(nOrb(k)+nDel(k))*(i-1+nFro(k)+nOcc(k))

do iSym=1,nSym
  !*********************************************************************
  ! Construct Wij(II)
  !*********************************************************************
  do iI=1,nOcc(iSym)
    E_i = EOcc(iOcc(iSym)+iI)
    do iJ=1,nFro(iSym)+nOcc(iSym)
      if (iJ <= nFro(iSym)) then
        E_j = EFro(iFro(iSym)+iJ)
      else
        E_j = EOcc(iOcc(iSym)+iJ-nFro(iSym))
      end if
      Work(iWDensActOcc(iI,iJ,iSym)) = Work(iWDensActOcc(iI,iJ,iSym))-Half*Work(iDensActOcc(iI,iJ,iSym))*(E_i+E_j)
    end do
  end do
  !*********************************************************************
  ! Construct Wab(II)
  !*********************************************************************
  do iA=1,nVir(iSym)
    E_a = EVir(iVir(iSym)+iA)
    do iB=1,nVir(iSym)+nDel(iSym)
      if (iB > nVir(iSym)) then
        E_b = EDel(iDel(iSym)+iB-nVir(iSym))
      else
        E_b = EVir(iVir(iSym)+iB)
      end if
      Work(iWDensVactVall(iA,iB,iSym)) = Work(iWDensVactVall(iA,iB,iSym))-Half*Work(iDensVactVall(iA,iB,iSym))*(E_a+E_b)
    end do
    !*******************************************************************
    ! Construct Wai(II) (The factor 2 in front of Pai is because Ppq is
    !                    already symmetrized here)
    !*******************************************************************
    do iI=1,nFro(iSym)+nOcc(iSym)
      if (iI <= nFro(iSym)) then
        E_i = EFro(iFro(iSym)+iI)
      else
        E_i = EOcc(iOcc(iSym)+iI-nFro(iSym))
      end if
      Work(iWDensVallOcc(iA,iI,iSym)) = Work(iWDensVallOcc(iA,iI,iSym))-Two*Work(iDensVallOcc(iA,iI,iSym))*(E_i)
    end do
  end do
end do

end subroutine Construct_WDensII
