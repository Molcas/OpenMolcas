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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      Subroutine  MemRg1(iAnga,nRys,MemPrm)
************************************************************************
*                                                                      *
* Modified to gradients 1991 R. Lindh, Dept. of Theoretical Chemistry, *
* University of Lund.                                                  *
************************************************************************
      Implicit Real*8 (a-h,o-z)
*
*     This routine will compute the memory requirement of Rysg1
*     Memory requirement is per primitive!
*
c#include "print.fh"
#include "itmax.fh"
      Integer iAnga(4)
*
*     Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
c     iRout = 13
c     iPrint = nPrint(iRout)
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      nRys = (la+lb+lc+ld+2 +1)/2
*
      MemPrm = 0
      nPAO = nElem(la)*nElem(lb)*nElem(lc)*nElem(ld)
*     1st order gradient of [ab|cd]
*     MemPrm = MemPrm + nPAO * 9
      nabMax = la+lb +1
      ncdMax = lc+ld +1
      nabcd = (nabMax+1)*(ncdMax+1)
      lB10=Max(Min(nabMax-1,1),0)
      lB01=Max(Min(ncdMax-1,1),0)
      lB00=Max(Min(Min(nabMax,ncdMax),1),0)
*     2D-Integrals
      n2D0=Max(nabcd,
     &         (la+2)*(lb+2)*(ncdMax+1),
     &         (la+2)*(lb+2)*(lc+2)*(ld+2))
      MemPrm = MemPrm + 3*nRys*n2D0
*     1st order gradient of the 2D-integrals
      n2D1=Max(nabcd,
     &         (la+2)*(lb+2)*(ncdMax+1),
     &         (la+1)*(lb+1)*(lc+1)*(ld+1)*3)
      MemPrm = MemPrm + 3*nRys*n2D1
*     Coefficients for recurrence relations
      MemPrm = MemPrm + 3*nRys + 3*nRys + 3*nRys*(lB10+lB01+lB00)
*     Roots
      MemPrm = MemPrm + nRys
*     The inverse of the arguments
      MemPrm = MemPrm + 1
*     Arguments
      MemPrm = MemPrm + 1
*     Expanded versions of Zeta, ZetInv, Eta, EtaInv,
*     P and Q
      MemPrm = MemPrm + 10
*     If (iPrint.ge.99) Then
*        Write (*,*) ' [ab|cd] 1st grad.   :', nPAO*9
*        Write (*,*) ' 2D-integrals        :', n2D0*3*nRys
*        Write (*,*) ' 2D-integrals (1st)  :', n2D1*3*nRys
*        Write (*,*) ' PAQP vector         :', 3*nRys
*        Write (*,*) ' QCPQ vector         :', 3*nRys
*        Write (*,*) ' B10 coefficients    :', nRys*3*lB10
*        Write (*,*) ' B00 coefficients    :', nRys*3*lB00
*        Write (*,*) ' B01 coefficients    :', nRys*3*lB01
*        Write (*,*) ' Roots               :', nRys
*        Write (*,*) ' Inverse arguments   :', 1
*        Write (*,*) ' Arguments           :', 1
*     End If
*
      Return
      End
