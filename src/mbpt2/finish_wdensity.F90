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
! Copyright (C) Jonas Bostrom                                          *
!***********************************************************************

subroutine Finish_WDensity()
! Author: Jonas Bostrom
!
! Purpose: Add the terms labeled [II] and [III] to the energy-weighted
!          MP2 Density

use MBPT2_Global, only: Density, EOcc, EVir, ip_Density, ip_WDensity, mAdDel, mAdFro, mAdOcc, mAdVir, WDensity
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iA, iB, iI, iJ, iP, iQ, iSym, iSym1, iSym2, iSymIJ, iSymPQ, nMaxOrb
real(kind=wp) :: Eps_a, Eps_b, Eps_i, Eps_j, Fac, xijpq, xipjq
real(kind=wp), allocatable :: Int1(:), IntC(:), Scr1(:)
#include "corbinf.fh"
! statement functions
integer(kind=iwp) :: i, j, k, iOccAOcc, iOccOcc, iVirVir, iVirOcc
iOccAOcc(i,j,k) = j-1+(nOrb(k)+nDel(k))*(nFro(k)+i-1)
iOccOcc(i,j,k) = j-1+(nOrb(k)+nDel(k))*(i-1)
iVirVir(i,j,k) = j-1+nFro(k)+nOcc(k)+(nOrb(k)+nDel(k))*(i-1+nFro(k)+nOcc(k))
iVirOcc(i,j,k) = j-1+nFro(k)+nOcc(k)+(nOrb(k)+nDel(k))*(i-1)

#ifdef _WARNING_WORKAROUND_
#include "compiler_features.h"
#if (( __GNUC__) && (GCC_VERSION < 70000))
iJ = 0
#endif
#endif

! Start with the II-terms
do iSym=1,nSym
  do iI=1,nOcc(iSym)
    Eps_i = EOcc(mAdOcc(iSym)+iI-1)
    do iJ=1,nFro(iSym)+nOcc(iSym)
      if (iJ <= nFro(iSym)) then
        Eps_j = EOcc(mAdFro(iSym)+iJ-1)
        Fac = Two
      else
        Eps_j = EOcc(mAdOcc(iSym)+iJ-nFro(iSym)-1)
        Fac = One
      end if
      ! Fac is because both core-occ and occ-core contributions should be
      ! added, they can both be added to IJ since we symmetrize W later and
      ! the term is symmetric
      WDensity(ip_WDensity(iSym)+iOccAOcc(iI,iJ,iSym)) = WDensity(ip_WDensity(iSym)+iOccAOcc(iI,iJ,iSym))- &
                                                         Fac*Density(ip_Density(iSym)+iOccAOcc(iI,iJ,iSym))*Half*(Eps_i+Eps_j)
    end do
  end do
  do iA=1,nExt(iSym)
    Eps_a = EVir(mAdVir(iSym)+iA-1)
    do iB=1,nExt(iSym)+nDel(iSym)
      if (iB > nExt(iSym)) then
        Eps_b = EVir(mAdDel(iSym)+iB-1-nExt(iSym))
      else
        Eps_b = EVir(mAdVir(iSym)+iB-1)
      end if
      WDensity(ip_WDensity(iSym)+iVirVir(iB,iA,iSym)) = WDensity(ip_WDensity(iSym)+iVirVir(iB,iA,iSym))- &
                                                        Density(ip_Density(iSym)+iVirVir(iB,iA,iSym))*Half*(Eps_a+Eps_b)
    end do
  end do

  do iI=1,nFro(iSym)+nOcc(iSym)
    if (iI <= nFro(iSym)) then
      Eps_i = EOcc(mAdFro(iSym)+iI-1)
    else
      Eps_i = EOcc(mAdOcc(iSym)+iI-nFro(iSym)-1)
    end if
    do iA=1,nExt(iSym)+nDel(iSym)
      WDensity(ip_WDensity(iSym)+iVirOcc(iI,iA,iSym)) = WDensity(ip_WDensity(iSym)+iVirOcc(iI,iA,iSym))- &
                                                        Density(ip_Density(iSym)+iVirOcc(iI,iA,iSym))*Two*(Eps_i)
    end do
  end do
end do
! And now the [III]-terms
nMaxOrb = 0
do iSym1=1,nSym
  do iSym2=1,nSym
    nMaxOrb = max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))*(nOrb(iSym2)+nDel(iSym2)))
  end do
end do
call mma_allocate(Int1,nMaxOrb,label='Int1')
call mma_allocate(IntC,nMaxOrb,label='IntC')
call mma_allocate(Scr1,nMaxOrb,label='Scr1')
do iSymIJ=1,nSym
  do iSymPQ=1,nSym
    do iI=1,nFro(iSymIJ)+nOcc(iSymIJ)
      do iJ=1,iI
        call Exch(iSymPQ,iSymIJ,iSymPQ,iSymIJ,iJ,iI,Int1,Scr1)
        call Coul(iSymPQ,iSymPQ,iSymIJ,iSymIJ,iJ,iI,IntC,Scr1)

        !write(u6,*) 'Finish'
        !write(u6,*) ' *  i,j = ',iI,iJ
        !call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSymPQ)+nDel(iSymPQ),nOrb(iSymPQ)+nDel(iSymPQ))
        !call RecPrt('IntC:','(8F10.6)',IntC,nOrb(iSymPQ)+nDel(iSymPQ),nOrb(iSymPQ)+nDel(iSymPQ))
        do iP=1,nOrb(iSymPQ)+nDel(iSymPQ)
          do iQ=1,nOrb(iSymPQ)+nDel(iSymPQ)
            xipjq = Int1(iP+(iQ-1)*(nOrb(iSymPQ)+nDel(iSymPQ)))
            xijpq = IntC(iP+(iQ-1)*(nOrb(iSymPQ)+nDel(iSymPQ)))
            WDensity(ip_WDensity(iSymIJ)+iOccOcc(iI,iJ,iSymIJ)) = &
              WDensity(ip_WDensity(iSymIJ)+iOccOcc(iI,iJ,iSymIJ))- &
              Density(ip_Density(iSymPQ)+iOccOcc(iP,iQ,iSymPQ))*(Two*xijpq-xipjq)
            if (iJ /= iI) then
              WDensity(ip_WDensity(iSymIJ)+iOccOcc(iJ,iI,iSymIJ)) = &
                WDensity(ip_WDensity(iSymIJ)+iOccOcc(iJ,iI,iSymIJ))- &
                Density(ip_Density(iSymPQ)+iOccOcc(iP,iQ,iSymPQ))*(Two*xijpq-xipjq)
            end if
          end do
        end do
      end do
    end do
  end do
end do
! Do the symmetrization of the energy weighted density and add the
! SCF-energyweighted density. Also changes sign on W since the formulas
! in the used article yields -W
do iSym=1,nSym
  do iP=1,nOrb(iSym)+nDel(iSym)
    do iQ=1,iP
      if (iQ /= iP) then
        WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -Half*(WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))+ &
                                                          WDensity(ip_WDensity(iSym)+iOccOcc(iQ,iP,iSym)))
        WDensity(ip_WDensity(iSym)+iOccOcc(iQ,iP,iSym)) = WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))
      else if (iQ <= nFro(iSym)) then
        WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))+ &
                                                          Two*EOcc(mAdFro(iSym)+iQ-1)
      else if (iQ <= nFro(iSym)+nOcc(iSym)) then
        WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))+ &
                                                          Two*EOcc(mAdOcc(iSym)+iQ-nFro(iSym)-1)
      else
        WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -WDensity(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))
      end if
    end do
  end do
end do

call mma_deallocate(Int1)
call mma_deallocate(IntC)
call mma_deallocate(Scr1)

return

end subroutine Finish_WDensity
