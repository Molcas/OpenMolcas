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

use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iA, iB, iI, iJ, iP, ipIntC, iQ, iSym, iSym1, iSym2, iSymIJ, iSymPQ, lint, nMaxOrb
real(kind=wp) :: Eps_a, Eps_b, Eps_i, Eps_j, Fac, xijpq, xipjq
#include "WrkSpc.fh"
#include "mp2grad.fh"
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
    Eps_i = Work(mAdOcc(iSym)+iI-1)
    do iJ=1,nFro(iSym)+nOcc(iSym)
      if (iJ <= nFro(iSym)) then
        Eps_j = Work(mAdFro(iSym)+iJ-1)
        Fac = Two
      else
        Eps_j = Work(mAdOcc(iSym)+iJ-nFro(iSym)-1)
        Fac = One
      end if
      ! Fac is because both core-occ and occ-core contributions should be
      ! added, they can both be added to IJ since we symmetrize W later and
      ! the term is symmetric
      Work(ip_WDensity(iSym)+iOccAOcc(iI,iJ,iSym)) = Work(ip_WDensity(iSym)+iOccAOcc(iI,iJ,iSym))-Fac*Work(ip_Density(iSym)+ &
                                                     iOccAOcc(iI,iJ,iSym))*Half*(Eps_i+Eps_j)
    end do
  end do
  do iA=1,nExt(iSym)
    Eps_a = Work(mAdVir(iSym)+iA-1)
    do iB=1,nExt(iSym)+nDel(iSym)
      if (iB > nExt(iSym)) then
        Eps_b = Work(mAdDel(iSym)+iB-1-nExt(iSym))
      else
        Eps_b = Work(mAdVir(iSym)+iB-1)
      end if
      Work(ip_WDensity(iSym)+iVirVir(iB,iA,iSym)) = Work(ip_WDensity(iSym)+iVirVir(iB,iA,iSym))-Work(ip_Density(iSym)+ &
                                                    iVirVir(iB,iA,iSym))*Half*(Eps_a+Eps_b)
    end do
  end do

  do iI=1,nFro(iSym)+nOcc(iSym)
    if (iI <= nFro(iSym)) then
      Eps_i = Work(mAdFro(iSym)+iI-1)
    else
      Eps_i = Work(mAdOcc(iSym)+iI-nFro(iSym)-1)
    end if
    do iA=1,nExt(iSym)+nDel(iSym)
      Work(ip_WDensity(iSym)+iVirOcc(iI,iA,iSym)) = Work(ip_WDensity(iSym)+iVirOcc(iI,iA,iSym))-Work(ip_Density(iSym)+ &
                                                    iVirOcc(iI,iA,iSym))*Two*(Eps_i)
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
lint = nMaxOrb
call GetMem('Int1','Allo','Real',ipInt1,lInt)
call GetMem('Int2','Allo','Real',ipInt2,lInt)
call GetMem('IntC','Allo','Real',ipIntC,lInt)
call GetMem('Scr1','Allo','Real',ipScr1,lInt)
do iSymIJ=1,nSym
  do iSymPQ=1,nSym
    do iI=1,nFro(iSymIJ)+nOcc(iSymIJ)
      do iJ=1,iI
        call Exch(iSymPQ,iSymIJ,iSymPQ,iSymIJ,iJ,iI,Work(ipInt1),Work(ipScr1))
        call Coul(iSymPQ,iSymPQ,iSymIJ,iSymIJ,iJ,iI,Work(ipIntC),Work(ipScr1))

        !write(u6,*) 'Finish'
        !write(u6,*) ' *  i,j = ',iI,iJ
        !call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSymPQ)+nDel(iSymPQ),nOrb(iSymPQ)+nDel(iSymPQ))
        !call RecPrt('IntC:','(8F10.6)',Work(ipIntC),nOrb(iSymPQ)+nDel(iSymPQ),nOrb(iSymPQ)+nDel(iSymPQ))
        do iP=1,nOrb(iSymPQ)+nDel(iSymPQ)
          do iQ=1,nOrb(iSymPQ)+nDel(iSymPQ)
            xipjq = Work(ipInt1+(iP-1)+(iQ-1)*(nOrb(iSymPQ)+nDel(iSymPQ)))
            xijpq = Work(ipIntC+(iP-1)+(iQ-1)*(nOrb(iSymPQ)+nDel(iSymPQ)))
            Work(ip_WDensity(iSymIJ)+iOccOcc(iI,iJ,iSymIJ)) = Work(ip_WDensity(iSymIJ)+iOccOcc(iI,iJ,iSymIJ))- &
                                                              Work(ip_Density(iSymPQ)+iOccOcc(iP,iQ,iSymPQ))*(Two*xijpq-xipjq)
            if (iJ /= iI) then
              Work(ip_WDensity(iSymIJ)+iOccOcc(iJ,iI,iSymIJ)) = Work(ip_WDensity(iSymIJ)+iOccOcc(iJ,iI,iSymIJ))- &
                                                                Work(ip_Density(iSymPQ)+iOccOcc(iP,iQ,iSymPQ))*(Two*xijpq-xipjq)
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
        Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -Half*(Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))+ &
                                                      Work(ip_WDensity(iSym)+iOccOcc(iQ,iP,iSym)))
        Work(ip_WDensity(iSym)+iOccOcc(iQ,iP,iSym)) = Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))
      else if (iQ <= nFro(iSym)) then
        Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))+Two*Work(mAdFro(iSym)+iQ-1)
      else if (iQ <= nFro(iSym)+nOcc(iSym)) then
        Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))+ &
                                                      Two*Work(mAdOcc(iSym)+iQ-nFro(iSym)-1)
      else
        Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))
      end if
    end do
  end do
end do

call GetMem('Int1','Free','Real',ipInt1,lInt)
call GetMem('Int2','Free','Real',ipInt2,lInt)
call GetMem('IntC','Free','Real',ipIntC,lInt)
call GetMem('Scr1','Free','Real',ipScr1,lInt)

return

end subroutine Finish_WDensity
