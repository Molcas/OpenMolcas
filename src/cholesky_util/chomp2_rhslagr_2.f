************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************

      SubRoutine ChoMP2_RHSlagr_2(EOcc,EVir,EFro,EDel,Xaibj,
     &                          LnPQRSprod,LiPQRSprod,iBatch,jBatch,
     &                          nOccLeftI,nOccLeftJ,
     &                          nOrbLeftI,nOrbLeftJ,
     &                          nFroLeftI,nFroLeftJ)
*     This will calculate the righthandside of the mp2lagrangian.
#include "implicit.fh"
      Real*8 EOcc(*), EVir(*),EFro(*),EDel(*), Xaibj(LnPQRSprod)
      Integer LiPQRSprod(8)
      Integer nOccLeftI(8),nOccLeftJ(8)
      Integer nOrbLeftI(8), nOrbLeftJ(8)
      Integer nFroLeftI(8),nFroLeftJ(8)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
*
      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      iDensVir(i,j,k) = ip_Density(k) +
     &                         nFro(k) + nOcc(k) +  j-1 +
     &                        (nOrb(k) + nDel(k))*
     &                        (i + nFro(k) + nOcc(k) -1)
      iDensOcc(i,j,k) = ip_Density(k)
     &                          + j-1 + (nOrb(k)+nDel(k))
     &                                      * (i-1)
      iMp2Lagr(i,j,k) = ip_Mp2Lagr(k) +
     &                          j-1 + (nOcc(k)+nFro(k))*(i-1)
      iDiaA(i,j,k) = ip_DiaA(k) + j-1 + (nOcc(k)+nFro(k))*(i-1)
      iFirstS(i,j)=iWork(ip_FirstS-1+nSym*(j-1)+i)
      LnBatOrb(i,j)=iWork(ip_LnBatOrb-1+nSym*(j-1)+i)
      LiPQprod(i,j,k)=iWork(ip_LiPQprod-1+
     &                        nSym*nSym*(k-1)+nSym*(j-1)+i)
      LnPQprod(i,j)=iWork(ip_LnPQprod-1+nSym*(j-1)+i)
*
      Do iSymIP = 1, nSym
         iSymAQ = iSymIP
         Do iSymI = 1,nSym
            iSymP = MulD2h(iSymI,iSymIP)
            LiOrb = Min(LnBatOrb(iSymI,iBatch),
     &                  nOccLeftI(iSymI))
            Do Li = 1, LiOrb
               If(nFroLeftI(iSymI) .gt. 0) Then
                  iI = Li
               Else
                  iI = iFirstS(iSymI,iBatch) + Li - 1
               End If
               Lii = LiPQprod(iSymI,iSymI,iBatch)
     &             + (nFro(iSymI) + nOcc(iSymI)
     &             +  nVir(iSymI) + nDel(iSymI))
     &             * (Li - 1)
     &             +  iI

*
               Do iSymA = 1, nSym
                  iSymAI = MulD2h(iSymI,iSymA)
                  iSymQ = MulD2h(iSymA,iSymAQ)
                  iSymIQ = MulD2h(iSymI,iSymQ)
                  LaOrb = LnBatOrb(iSymA,jBatch)
     &                  - nOccLeftJ(iSymA)
                  Do La = 1, LaOrb
                     If(nOccLeftJ(iSymA).gt.0) Then
                        iA = La
                     Else
                        iA = iFirstS(iSymA,jBatch) - nFro(iSymA)
     &                     - nOcc(iSymA) + La-1
                     End If
*     Calculating the diagonal of A: 3(ia|ia) - (ii|aa)
*     (Only needed once for each iA hence iP = 1)
*
                     If((iSymI .eq. iSymA).and.(iSymIP.eq.1)) Then
                        Lii = LiPQprod(iSymI,iSymI,iBatch)
     &                       + (nFro(iSymI) + nOcc(iSymI)
     &                       +  nVir(iSymI) + nDel(iSymI))
     &                       * (Li - 1)
     &                       +  iI
                        Laa = LiPQprod(iSymA,iSymA,jBatch)
     &                       + (nFro(iSymA) + nOcc(iSymA)
     &                       +  nVir(iSymA) + nDel(iSymA))
     &                       * (La + nOccLeftJ(iSymA) - 1)
     &                       +  iA + nFro(iSymA) + nOcc(iSymA)
                        Lai = LiPQprod(iSymA,iSymI,iBatch)
     &                       + (nFro(iSymA) + nOcc(iSymA)
     &                       +  nVir(iSymA) + nDel(iSymA))
     &                       * (Li-1)
     &                       +  iA + nFro(iSymA) + nOcc(iSymA)
                        Lia = LiPQprod(iSymI,iSymA,jBatch)
     &                       + (nFro(iSymA) + nOcc(iSymA)
     &                       +  nVir(iSymA) + nDel(iSymA))
     &                       * (La + nOccLeftJ(iSymA) - 1)
     &                       +  iI
                        If(iBatch.eq.jBatch) Then
                           ip_iiaa = LiPQRSprod(iSymAI)
     &                             + iTri(Laa,Lii)
                           ip_aiai = LiPQRSprod(iSymAI)
     &                             + iTri(Lia,Lai)
                        Else
                           ip_iiaa = LiPQRSprod(iSymAI)
     &                             + LnPQprod(iSymAI,iBatch)
     &                             * (Laa-1) + Lii
                           ip_aiai = LiPQRSprod(iSymAI)
     &                             + LnPQprod(iSymAI,iBatch)
     &                             * (Lia-1) + Lai
                        End If
                        DiagInt = 3.0d0 * Xaibj(ip_aiai)
     &                                  - Xaibj(ip_iiaa)
                        If(iI .le. nFro(iSymI)) Then
                           E_Occ = EFro(iFro(iSymI) + iI)
                        Else
                           E_Occ = EOcc(iOcc(iSymI)
     &                           + iI - nFro(iSymI))
                        End If
                        If(iA .gt. nVir(iSymA)) Then
                           E_Vir = EDel(iDel(iSymA)
     &                           + iA - nVir(iSymA))
                        Else
                           E_Vir = EVir(iVir(iSymA) + iA)
                        End If

                        DiagEn = E_Vir - E_Occ
                        Work(iDiaA(iA,iI,iSymI)) =
     &                       Work(iDiaA(iA,iI,iSymI))
     &                       + 1.0D0/(DiagInt + DiagEn)
*-----------------------------------------------------------------------
*                              Write(6,*) 'IA', iI, iA
*                              Write(6,*) 'Symm', iSymI, iSymA
*                              Write(6,*) 'iiaa', Xaibj(ip_iiaa)
*                              Write(6,*) 'iaia', Xaibj(ip_aiai)
*                              Write(6,*) 'DiagInt', DiagInt
*                              Write(6,*) 'DiagEn', DiagEn
*-----------------------------------------------------------------------
                     End If
                     Do iP = 1, nOrb(iSymP) + nDel(iSymP)
                        Lip = LiPQprod(iSymP,iSymI,iBatch)
     &                       + (nFro(iSymP) + nOcc(iSymP)
     &                       +  nVir(iSymP) + nDel(iSymP))
     &                       * (Li - 1)
     &                       +  iP
                        Lap = LiPQProd(iSymP,iSymA,jBatch)
     &                       + (nFro(iSymP) + nOcc(iSymP)
     &                       +  nVir(iSymP) + nDel(iSymP))
     &                       * (La+nOccLeftJ(iSymA)-1)
     &                       +  iP
                        If((iP .le. nFro(iSymP) + nOcc(iSymP)) .and.
     &                      (iSymI .eq. iSymP)) Then
                           Do iK = 1, nFro(iSymQ) + nOcc(iSymQ)
                              iJ = iP
                              Lak = LiPQProd(iSymQ,iSymA,jBatch)
     &                            + (nFro(iSymQ) + nOcc(iSymQ)
     &                            +  nVir(iSymQ) + nDel(iSymQ))
     &                            * (La+nOccLeftJ(iSymA)-1)
     &                            +  iK
                              Lik = LiPQProd(iSymQ,iSymI,iBatch)
     &                             + (nFro(iSymQ) + nOcc(iSymQ)
     &                             +  nVir(iSymQ) + nDel(iSymQ))
     &                             * (Li-1)
     &                             +  iK
*
                              If(iBatch.eq.jBatch) Then
                                 ip_ijak = LiPQRSprod(iSymIP)
     &                                + iTri(Lip,Lak)
                                 ip_ikaj = LiPQRSprod(iSymIQ)
     &                             + iTri(Lik,Lap)
                              Else
                                 ip_ijak = LiPQRSprod(iSymIP)
     &                                + LnPQprod(iSymIP,iBatch)
     &                                * (Lak-1) + Lip
                                 ip_ikaj = LiPQRSprod(iSymIQ)
     &                                + LnPQprod(iSymIQ,iBatch)
     &                                * (Lap-1) + Lik
                              End If
                              A = 2.0d0 * Xaibj(ip_ijak)
     &                          - 1.0d0 * Xaibj(ip_ikaj)
                              Work(iMp2Lagr(iA,iK,iSymQ)) =
     &                             Work(iMp2Lagr(iA,iK,iSymQ))
     &                             + Work(iDensOcc(iJ,iI,iSymI))
     &                             * A
*--------------------- Debug Comments -----------------------
*                              If(iP .le. nFro(iSymP) + nOcc(iSymP)) Then
*                                 Write(6,*) 'AIJK',iA,iI,iJ,iK
*                                 Write(6,*) 'aijk', Xaibj(ip_ijak)
*                                 Write(6,*) 'akji', Xaibj(ip_ikaj)
*                                 Write(6,*) 'A',A
*                                 Write(6,*) 'Density',Work(iDensOcc(iJ,
*     &                                                     iI,iSymI))
*                              End If
*------------------------------------------------------------
                           End Do
                        Else If((iP .gt. nFro(iSymP) + nOcc(iSymP)).and.
     &                          (iSymA .eq. iSymQ)) Then
                           Do iB = 1, nVir(iSymQ) + nDel(iSymQ)
                              iC = iP - nFro(iSymP) - nOcc(iSymP)
                              Lab = LiPQProd(iSymQ,iSymA,jBatch)
     &                            + (nFro(iSymQ) + nOcc(iSymQ)
     &                            +  nVir(iSymQ) + nDel(iSymQ))
     &                            * (La+nOccLeftJ(iSymA)-1)
     &                            +  iB + nFro(iSymQ) + nOcc(iSymQ)
                              Lib = LiPQProd(iSymQ,iSymI,iBatch)
     &                             + (nFro(iSymQ) + nOcc(iSymQ)
     &                             +  nVir(iSymQ) + nDel(iSymQ))
     &                             * (Li-1)
     &                             + iB + nFro(iSymQ) + nOcc(iSymQ)
                              If(iBatch.eq.jBatch) Then
                                 ip_icab = LiPQRSprod(iSymIP)
     &                                + iTri(Lip,Lab)
                                 ip_ibac = LiPQRSprod(iSymIQ)
     &                             + iTri(Lib,Lap)
                              Else
                                 ip_icab = LiPQRSprod(iSymIP)
     &                                + LnPQprod(iSymIP,iBatch)
     &                                * (Lab-1) + Lip
                                 ip_ibac = LiPQRSprod(iSymIQ)
     &                                + LnPQprod(iSymIQ,iBatch)
     &                                * (Lap-1) + Lib
                              End If
                              A = 2.0d0 * Xaibj(ip_icab)
     &                          - 1.0d0 * Xaibj(ip_ibac)
                              Work(iMp2Lagr(iC,iI,iSymI)) =
     &                             Work(iMp2Lagr(iC,iI,iSymI))
     &                             + Work(iDensVir(iB,iA,iSymA))
     &                             * A
*----------- Debug Comments ------------------------------------------
*                        Write(6,*) 'AIBC',iA,iI,iB,iC
*                        Write(6,*) 'Symm',iSymA,iSymI,iSymQ,iSymP
*                        Write(6,*) 'icab', Xaibj(ip_icab)
*                        Write(6,*) 'ibac', Xaibj(ip_ibac)
*                        Write(6,*) 'adress',iMp2Lagr(iC,iI,iSymI)-
*     &                             ip_Mp2Lagr(1)
*                        Write(6,*) 'DensAdress',iDensVir(iB,iA,iSymA)-
*     &                                          ip_Density(1)
*                        Write(6,*) 'Dens',Work(iDensVir(iB,iA,iSymA))
*                        Write(6,*) 'A',A
*                        Write(6,*) 'Bidrag',
*     &                             A*Work(iDensVir(iB,iA,iSymA))
*----------------------------------------------------------------------
                           End Do
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(nOrbLeftI)
         Call Unused_integer_array(nOrbLeftJ)
         Call Unused_integer_array(nFroLeftJ)
      End If
      End
