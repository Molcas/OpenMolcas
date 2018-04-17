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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Subroutine  MemRys(iAnga,MemPrm)
      Implicit Real*8 (a-h,o-z)
*
*     This routine will compute the memory requirement of RYS
*     Memory requirement is per primitive!
*
#include "itmax.fh"
#include "print.fh"
#include "FMM.fh"
cgh - stuff for short range integrals
#include "srint.fh"
      Integer iAnga(4)
*
*     Statement function for canonical index, etc.
*
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 13
      iPrint = nPrint(iRout)
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      nRys = (la+lb+lc+ld+2)/2
      labMin=nabSz(Max(la,lb)-1)+1
      labMax=nabSz(la+lb)
      lcdMin=nabSz(Max(lc,ld)-1)+1
      lcdMax=nabSz(lc+ld)
      labcd = (labMax-labMin+1)*(lcdMax-lcdMin+1)
      If (iPrint.ge.99) Then
         Write (6,*) ' labMin=',labMin
         Write (6,*) ' labMax=',labMax
         Write (6,*) ' lcdMin=',lcdMin
         Write (6,*) ' lcdMax=',lcdMax
      End If
      MemPrm = 0
*     [a0|c0]
      MemPrm = MemPrm + labcd
*                                                                      *
************************************************************************
*                                                                      *
*     For FMM, we only want short-range integrals, using twice the memory
*     to store full and long-range components (which are subtracted)
*    -same for MOLPRO shortrange
*
      If (FMM_shortrange.or.shortrange) MemPrm = MemPrm + labcd
*                                                                      *
************************************************************************
*                                                                      *
      nabMax = la+lb
*     nabMin = Max(la,lb)
      ncdMax = lc+ld
*     ncdMin = Max(lc,ld)
      nabcd = (nabMax+1)*(ncdMax+1)
      lB10=Max(Min(nabMax-1,1),0)
      lB01=Max(Min(ncdMax-1,1),0)
      lB00=Max(Min(Min(nabMax,ncdMax),1),0)
*     Normalization
      MemPrm = MemPrm + 1
*     2D-Integrals
      MemPrm = MemPrm + nabcd*3*nRys
*     Coefficients for recurrence relations
      MemPrm = MemPrm + 3*nRys + 3*nRys + 3*nRys*(lB10+lB01+lB00)
*     Roots
      MemPrm = MemPrm + nRys
*     The inverse of the arguments
      MemPrm = MemPrm + 1
*     Arguments
      MemPrm = MemPrm + 1
*     Expanded versions of Zeta, ZetInv, Eta, EtaInv,
*     rKapab, rKapcd, P and Q
      MemPrm = MemPrm + 12
      If (iPrint.ge.99) Then
         Write (6,*) ' [e0|f0] integrals   :', labcd
         Write (6,*) ' Normalization factor:', 1
         Write (6,*) ' 2D-integrals        :', nabcd*3*nRys
         Write (6,*) ' PAQP vector         :', 3*nRys
         Write (6,*) ' QCPQ vector         :', 3*nRys
         Write (6,*) ' B10 coefficients    :', nRys*3*lB10
         Write (6,*) ' B00 coefficients    :', nRys*3*lB00
         Write (6,*) ' B01 coefficients    :', nRys*3*lB01
         Write (6,*) ' Roots               :', nRys
         Write (6,*) ' Inverse arguments   :', 1
         Write (6,*) ' Arguments           :', 1
      End If
*
      Return
      End
