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
* Copyright (C) 2008, Jonas Bostrom                                    *
************************************************************************

      Subroutine Construct_WDensII(EOcc,EVir,EFro,EDel)
*
*     Jonas Bostrom, October 2008
*
*     Purpose: Construct the piece of the energy-weighted density
*              usually labeled II.
*
#include "implicit.fh"
      Real*8 EOcc(*), EVir(*), EFro(*), EDel(*)
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
*
      Character*20 ThisNm
      Character*25 SecNam
      Parameter (SecNam = 'Construct_WDensII',
     &           ThisNm = 'Construct_WDensII')
*
      iDensActOcc(i,j,k) = ip_Density(k)
     &                   +  j-1
     &                   + (nOrb(k) + nDel(k))
     &                   * (i + nFro(k) - 1)
      iWDensActOcc(i,j,k) = ip_WDensity(k)
     &                    +  j-1
     &                    + (nOrb(k) + nDel(k))
     &                    * (i-1 + nFro(k))
      iDensVactVall(i,j,k) = ip_Density(k)
     &                     + j-1 + nFro(k) + nOcc(k)
     &                     + (nOrb(k) + nDel(k))
     &                     * (i-1 + nFro(k) + nOcc(k))
      iWDensVactVall(i,j,k) = ip_WDensity(k)
     &                      + j-1 + nFro(k) + nOcc(k)
     &                      + (nOrb(k) + nDel(k))
     &                      * (i-1 + nFro(k) + nOcc(k))
      iDensVallOcc(i,j,k) = ip_Density(k)
     &                    +  j-1
     &                    + (nOrb(k) + nDel(k))
     &                    * (i-1 + nFro(k) + nOcc(k))
      iWDensVallOcc(i,j,k) = ip_WDensity(k)
     &                    +  j-1
     &                    + (nOrb(k) + nDel(k))
     &                    * (i-1 + nFro(k) + nOcc(k))
*
      Do iSym = 1, nSym
********************************************************************
*     Construct Wij(II)
********************************************************************
         Do iI = 1, nOcc(iSym)
            E_i = EOcc(iOcc(iSym) + iI)
            Do iJ = 1, nFro(iSym) + nOcc(iSym)
               If(iJ .le. nFro(iSym)) Then
                  E_j = EFro(iFro(iSym) + iJ)
               Else
                  E_j = EOcc(iOcc(iSym) + iJ - nFro(iSym))
               End If
               Work(iWDensActOcc(iI,iJ,iSym)) =
     &              Work(iWDensActOcc(iI,iJ,iSym))
     &           -  0.5d0 * Work(iDensActOcc(iI,iJ,iSym))
     &           * (E_i + E_j)
            End Do
         End Do
*********************************************************************
*     Construct Wab(II)
*********************************************************************
         Do iA = 1, nVir(iSym)
            E_a = EVir(iVir(iSym) + iA)
            Do iB = 1, nVir(iSym) + nDel(iSym)
               If(iB .gt. nVir(iSym)) Then
                  E_b = EDel(iDel(iSym) + iB - nVir(iSym))
               Else
                  E_b = EVir(iVir(iSym) + iB)
               End If
               Work(iWDensVactVall(iA,iB,iSym)) =
     &              Work(iWDensVactVall(iA,iB,iSym))
     &           -  0.5d0 * Work(iDensVactVall(iA,iB,iSym))
     &           * (E_a + E_b)
            End Do
************************************************************************
*     Construct Wai(II) (The factor 2 in front of Pai is because Ppq is
*                        already symmetrized here)
************************************************************************
            Do iI = 1, nFro(iSym) + nOcc(iSym)
               If(iI .le. nFro(iSym)) Then
                  E_i = EFro(iFro(iSym) + iI)
               Else
                  E_i = EOcc(iOcc(iSym) + iI - nFro(iSym))
               End If
               Work(iWDensVallOcc(iA,iI,iSym)) =
     &              Work(iWDensVallOcc(iA,iI,iSym))
     &           - 2.0d0 * Work(iDensVallOcc(iA,iI,iSym))
     &           * (E_i)
            End Do
         End Do
      End Do
      End
