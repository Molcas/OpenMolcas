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

      SubRoutine Construct_WDensIII(Xaibj,LnPQRSprod,LiPQRSprod,iBatch,
     &                              jBatch,nOccLeftI,nOccLeftJ)
*
*     Jonas Bostrom, October 2008
*
*     Purpose: Construct the piece of the energy-weighted density
*              usually labeled III.
*
#include "implicit.fh"
      Real*8 Xaibj(LnPQRSprod)
      Integer LiPQRSprod(8)
      Integer nOccLeftI(8),nOccLeftJ(8)
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
*
      Character*20 ThisNm
      Character*25 SecNam
      Parameter (SecNam = 'Construct_WDensIII',
     &           ThisNm = 'Construct_WDensIII')
*
      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      iDensAllAll(i,j,k) = ip_Density(k)
     &                   +  j-1
     &                   + (nOrb(k) + nDel(k))
     &                   * (i-1)
      iWDensOccOcc(i,j,k) = ip_WDensity(k)
     &                    +  j-1
     &                    + (nOrb(k) + nDel(k))
     &                    * (i-1)
      iFirstS(i,j)=iWork(ip_FirstS-1+nSym*(j-1)+i)
      LnBatOrb(i,j)=iWork(ip_LnBatOrb-1+nSym*(j-1)+i)
      LiPQprod(i,j,k) = iWork(ip_LiPQprod-1+
     &                        nSym*nSym*(k-1)+nSym*(j-1)+i)
      LnPQprod(i,j)=iWork(ip_LnPQprod-1+nSym*(j-1)+i)
*
      iSymPQ = 1
      iSymIJ = 1
      Do iSymJ = 1, nSym
         LjOrb = Min(LnBatOrb(iSymJ,jBatch),
     &                  nOccLeftJ(iSymJ))
         Do Lj = 1, LjOrb
            iJ = iFirstS(iSymJ,jBatch) + Lj-1
            Do iSymQ = 1,nSym
               Do iQ = 1, nOrb(iSymQ) + nDel(iSymQ)
                  iSymP = MulD2h(iSymPQ,iSymQ)
                  iSymQJ = MulD2h(iSymQ,iSymJ)
                  iSymPI = iSymQJ
                  iSymP = MulD2h(iSymPQ,iSymQ)
                  Ljq = LiPQProd(iSymQ,iSymJ,jBatch)
     &                + (nFro(iSymQ) + nOcc(iSymQ)
     &                +  nVir(iSymQ) + nDel(iSymQ))
     &                * (Lj-1)
     &                + iQ
                  LpOrb = LnBatOrb(iSymP,iBatch)
                  Do Lp = 1, LpOrb
                     iP = iFirstS(iSymP,iBatch) + Lp-1
                     iSymI = MulD2h(iSymIJ,iSymJ)
                     Lpq = LiPQProd(iSymQ,iSymP,iBatch)
     &                   + (nFro(iSymQ) + nOcc(iSymQ)
     &                   +  nVir(iSymQ) + nDel(iSymQ))
     &                   * (Lp-1)
     &                   + iQ
                     Do iI = 1, nFro(iSymI) + nOcc(iSymI)
                        Lpi = LiPQProd(iSymI,iSymP,iBatch)
     &                      + (nFro(iSymI) + nOcc(iSymI)
     &                      +  nVir(iSymI) + nDel(iSymI))
     &                      * (Lp-1)
     &                      + iI
                        Lji = LiPQProd(iSymI,iSymJ,jBatch)
     &                      + (nFro(iSymI) + nOcc(iSymI)
     &                      +  nVir(iSymI) + nDel(iSymI))
     &                      * (Lj-1)
     &                      + iI
*
                        If(iBatch.eq.jBatch) Then
                           ip_pqij = LiPQRSprod(iSymPQ)
     &                             + iTri(Lji,Lpq)
                           ip_piqj = LiPQRSprod(iSymPI)
     &                             + iTri(Ljq,Lpi)
                        Else
                           ip_pqij = LiPQRSprod(iSymPQ)
     &                             + LnPQprod(iSymPQ,iBatch)
     &                             * (Lji-1) + Lpq
                           ip_piqj = LiPQRSprod(iSymPI)
     &                             + LnPQprod(iSymPI,iBatch)
     &                             * (Ljq-1) + Lpi
                        End If
*
                        X = 2.0d0 * Xaibj(ip_pqij) - Xaibj(ip_piqj)
                        Work(iWDensOccOcc(iI,iJ,iSymJ)) =
     &                       Work(iWDensOccOcc(iI,iJ,iSymJ))
     &                     - X * Work(iDensAllAll(iP,iQ,iSymQ))
*----------- Debug Comments -------------------------------------------
*                        Write(6,*) 'IJPQ', iI, iJ, iP, iQ
*                        Write(6,*) 'Symm',iSymI, iSymJ, iSymP, iSymQ
*                        Write(6,*) 'pqij', Xaibj(ip_pqij)
*                        Write(6,*) 'piqj', Xaibj(ip_piqj)
*                        Write(6,*) 'Dens',Work(iDensAllAll(iP,iQ,iSymQ))
*-----------------------------------------------------------------------

                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nOccLeftI)
      End
