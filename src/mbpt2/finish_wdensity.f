************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Jonas Bostrom                                          *
************************************************************************
      subroutine Finish_WDensity()
*********************************************
*
*     Author: Jonas Bostrom
*
*     Purpose: Add the terms labeled [II] and [III] to the energy-weighted
*              MP2 Density
*
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "files_mbpt2.fh"
#include "trafo.fh"
#include "mp2grad.fh"
#include "corbinf.fh"
      iOccAOcc(i,j,k) =   j-1 + (nOrb(k)+nDel(k)) *
     &                              (nFro(k) + i-1)
      iOccOcc(i,j,k) = j-1 + (nOrb(k) + nDel(k)) * (i-1)
      iVirVir(i,j,k) = j-1 + nFro(k) + nOcc(k) +
     &                      (nOrb(k) + nDel(k))*(i-1+nFro(k) + nOcc(k))
      iVirOcc(i,j,k) = j-1 + nFro(k) + nOcc(k) + (nOrb(k) + nDel(k))
     &                                *(i-1)
*
*     Get rid of silly compiler warning.
      iJ = 0
*     Start with the II-terms
      Do iSym = 1, nSym
         Do iI = 1, nOcc(iSym)
            Eps_i = Work(mAdOcc(iSym)+iI-1)
            Do iJ = 1, nFro(iSym) + nOcc(iSym)
               If(iJ.le.nFro(iSym)) Then
                  Eps_j = Work(mAdFro(iSym)+iJ-1)
                  Fac = 2.0d0
               Else
                  Eps_j = Work(mAdOcc(iSym)+iJ-nFro(iSym)-1)
                  Fac = 1.0d0
               End If
*     Fac is because both core-occ and occ-core contributions should be
*     added, they can both be added to IJ since we symmetrize W later and
*     the term is symmetric
               Work(ip_WDensity(iSym) + iOccAOcc(iI,iJ,iSym)) =
     &              Work(ip_WDensity(iSym) + iOccAOcc(iI,iJ,iSym)) -
     &              Fac*Work(ip_Density(iSym) + iOccAOcc(iI,iJ,iSym))*
     &                  0.5d0*(Eps_i + Eps_j)
             End Do
         End Do
         Do iA = 1, nExt(iSym)
            Eps_a = Work(mAdVir(iSym) + iA-1)
            Do iB = 1, nExt(iSym) + nDel(iSym)
               If(iB.gt.nExt(iSym)) Then
                  Eps_b = Work(mAdDel(iSym) + iB-1-nExt(iSym))
               Else
                  Eps_b = Work(mAdVir(iSym) + iB-1)
               End If
               Work(ip_WDensity(iSym) + iVirVir(iB,iA,iSym)) =
     &              Work(ip_WDensity(iSym) + iVirVir(iB,iA,iSym)) -
     &              Work(ip_Density(iSym) + iVirVir(iB,iA,iSym))*
     &              0.5d0*(Eps_a + Eps_b)
            End Do
         End Do

         Do iI = 1, nFro(iSym) + nOcc(iSym)
            If(iI.le.nFro(iSym)) Then
               Eps_i = Work(mAdFro(iSym)+iI-1)
            Else
               Eps_i = Work(mAdOcc(iSym)+iI-nFro(iSym)-1)
            End If
            Do iA = 1, nExt(iSym) + nDel(iSym)
               Work(ip_WDensity(iSym) + iVirOcc(iI,iA,iSym)) =
     &              Work(ip_WDensity(iSym) + iVirOcc(iI,iA,iSym)) -
     &              Work(ip_Density(iSym) + iVirOcc(iI,iA,iSym))*
     &              2.0d0*(Eps_i)
            End Do
         End Do
      End Do
*     And now the [III]-terms
      nMaxOrb=0
      Do iSym1 = 1,nSym
         Do iSym2 = 1,nSym
            nMaxOrb = Max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))
     &              *(nOrb(iSym2)+nDel(iSym2)))
         End Do
      End Do
      lint = nMaxOrb
      Call GetMem('Int1','Allo','Real',ipInt1,lInt)
      Call GetMem('Int2','Allo','Real',ipInt2,lInt)
      Call GetMem('IntC','Allo','Real',ipIntC,lInt)
      Call GetMem('Scr1','Allo','Real',ipScr1,lInt)
      Do iSymIJ = 1, nSym
         Do iSymPQ = 1, nSym
            Do iI = 1, nFro(iSymIJ) + nOcc(iSymIJ)
               Do iJ = 1, iI
                  Call Exch(iSymPQ,iSymIJ,iSymPQ,iSymIJ,
     &                 iJ,
     &                 iI,
     &                 Work(ipInt1), Work(ipScr1))
                  Call Coul(iSymPQ,iSymPQ,iSymIJ,iSymIJ,
     &                 iJ,
     &                 iI,
     &                 Work(ipIntC), Work(ipScr1))

*                  Write(6,*) 'Finish'
*                  Write(6,*) ' *  i,j = ',iI, iJ
*               Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),
*     &              nOrb(iSymPQ)+nDel(iSymPQ),
*     &              nOrb(iSymPQ)+nDel(iSymPQ))
*               Call RecPrt('IntC:','(8F10.6)',Work(ipIntC),
*     &              nOrb(iSymPQ)+nDel(iSymPQ),
*     &              nOrb(iSymPQ)+nDel(iSymPQ))
                  Do iP = 1, nOrb(iSymPQ) + nDel(iSymPQ)
                     Do iQ = 1, nOrb(iSymPQ)+nDel(iSymPQ)
                        xipjq = Work(ipInt1 + (iP-1) +
     &                         (iQ-1)*(nOrb(iSymPQ) + nDel(iSymPQ)))
                        xijpq = Work(ipIntC + (iP-1) +
     &                         (iQ-1)*(nOrb(iSymPQ) + nDel(iSymPQ)))
                        Work(ip_WDensity(iSymIJ)+iOccOcc(iI,iJ,iSymIJ))=
     &                  Work(ip_WDensity(iSymIJ)+iOccOcc(iI,iJ,iSymIJ))-
     &                  Work(ip_Density(iSymPQ)+iOccOcc(iP,iQ,iSymPQ))*
     &                   (2.0d0*xijpq - xipjq)
                        If(iJ.ne.iI) Then
                           Work(ip_WDensity(iSymIJ)
     &                         +iOccOcc(iJ,iI,iSymIJ))=
     &                          Work(ip_WDensity(iSymIJ)
     &                          +iOccOcc(iJ,iI,iSymIJ))-
     &                          Work(ip_Density(iSymPQ)
     &                          +iOccOcc(iP,iQ,iSymPQ))*
     &                          (2.0d0*xijpq - xipjq)
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
*     Do the symmetrization of the energy weighted density and add the
*     SCF-energyweighted density. Also changes sign on W since the formulas
*     in the used article yields -W
      Do iSym = 1, nSym
         Do iP=1, nOrb(iSym) + nDel(iSym)
            Do iQ =1,iP
               If(iQ.ne.iP) Then
                  Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -
     &             (Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) +
     &              Work(ip_WDensity(iSym)+iOccOcc(iQ,iP,iSym)))/2.0d0
                  Work(ip_WDensity(iSym)+iOccOcc(iQ,iP,iSym)) =
     &                Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))
               Else If(iQ.le.nFro(iSym)) Then
                  Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -
     &              Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) +
     &              2.0d0*Work(mAdFro(iSym)+iQ-1)
               Else If(iQ.le.nFro(iSym)+nOcc(iSym)) Then
                  Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -
     &              Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) +
     &              2.0d0*Work(mAdOcc(iSym)+iQ-nFro(iSym)-1)
               Else
                  Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym)) = -
     &              Work(ip_WDensity(iSym)+iOccOcc(iP,iQ,iSym))
               End If
            End Do
         End Do
      End Do

*
      Call GetMem('Int1','Free','Real',ipInt1,lInt)
      Call GetMem('Int2','Free','Real',ipInt2,lInt)
      Call GetMem('IntC','Free','Real',ipIntC,lInt)
      Call GetMem('Scr1','Free','Real',ipScr1,lInt)

      Return
      End
