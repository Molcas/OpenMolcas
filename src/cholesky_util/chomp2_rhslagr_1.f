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

      SubRoutine ChoMP2_rhslagr_1(EOcc,EVir,EFro,EDel,Xaibj,
     &                            LnPQRSprod,LiPQRSprod,iBatch,jBatch,
     &                            nOccLeftI,nOccLeftJ,
     &                            nOrbLeftI,nOrbLeftJ,
     &                            nFroLeftI,nFroLeftJ)
*     This will calculate the righthandside of the mp2lagrangian.
*
#include "implicit.fh"
      Real*8 EOcc(*), EVir(*),EFro(*),EDel(*), Xaibj(LnPQRSprod)
      Integer LiPQRSprod(8)
      Integer nOccLeftI(8),nOccLeftJ(8)
      Integer nOrbLeftI(8),nOrbLeftJ(8)
      Integer nFroLeftI(8),nFroLeftJ(8)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "WrkSpc.fh"
*
      Character*10  ThisNm
      Character*17 SecNam
      Parameter (SecNam = 'ChoMP2_rhslagr_1', ThisNm = 'rhslagr_1')
*
      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      iDensVir(i,j,k) = ip_Density(k) +
     &                         nFro(k) + nOcc(k) +  j-1 +
     &                        (nOrb(k) + nDel(k))*
     &                        (i + nFro(k) + nOcc(k) -1)
      iWDensVir(i,j,k) = ip_WDensity(k) +
     &                         nFro(k) + nOcc(k) + j-1
     &                      + (nOrb(k) + nDel(k))
     &                      * (i + nFro(k) + nOcc(k) - 1)
      iDensVirFro(i,j,k) = ip_Density(k) +
     &                         nFro(k) + nOcc(k) + nVir(k) + j-1 +
     &                        (nOrb(k) + nDel(k))*
     &                        (i + nFro(k) + nOcc(k) -1)
      iWDensVirFro(i,j,k) = ip_WDensity(k) +
     &                         nFro(k) + nOcc(k) + nVir(k) + j-1 +
     &                        (nOrb(k) + nDel(k))*
     &                        (i + nFro(k) + nOcc(k) -1)
      iDensFroVir(i,j,k) = ip_Density(k) +
     &                         nFro(k) + nOcc(k) + j-1 +
     &                        (nOrb(k) + nDel(k))
     &                      * (i + nFro(k) + nOcc(k) + nVir(k) - 1)
      iWDensOcc(i,j,k) = ip_WDensity(k)
     &                    +  j + nFro(k) -1
     &                    + (nOrb(k) + nDel(k))
     &                    * (i + nFro(k) -1)
      iWDensOccFro(i,j,k) = ip_WDensity(k)
     &                    +  j-1
     &                    + (nOrb(k) + nDel(k))
     &                    * (i + nFro(k) -1)

      iDensOcc(i,j,k) = ip_Density(k) +
     &                        nFro(k) + j-1 + (nOrb(k)+nDel(k))
     &                                      * (i + nFro(k)-1)
      iDensFroOcc(i,j,k) = ip_Density(k) +
     &                            j-1 + nFro(k) + (nOrb(k)+nDel(k))
     &                                * (i-1)
      iDensOccFro(i,j,k) = ip_Density(k) +
     &                            j-1 + (nOrb(k)+nDel(k))
     &                                * (i + nFro(k) -1)
      iMp2Lagr(i,j,k) = ip_Mp2Lagr(k) +
     &                          j-1 + (nOcc(k)+nFro(k))*(i-1)
      iWDensVactOcc(i,j,k) = ip_WDensity(k) +
     &                           j-1
     &                        + (nOrb(k) + nDel(k))
     &                        * (i + nFro(k) + nOcc(k) - 1)
      iFirstS(i,j)=iWork(ip_FirstS-1+nSym*(j-1)+i)
      LnBatOrb(i,j)=iWork(ip_LnBatOrb-1+nSym*(j-1)+i)
      LiPQprod(i,j,k) = iWork(ip_LiPQprod-1+
     &                        nSym*nSym*(k-1)+nSym*(j-1)+i)
      LnPQprod(i,j)=iWork(ip_LnPQprod-1+nSym*(j-1)+i)
*
      Do iSymBJ = 1,nSym
         iSymAI = iSymBJ
********************************************************************
****   Common code for Pab, Wab, Lagr(3) and PaB
********************************************************************
         Do iSymJ = 1,nSym
            iSymB = MulD2h(iSymJ,iSymBJ)
            LjOrb = Min(LnBatOrb(iSymJ,jBatch)-nFroLeftJ(iSymJ),
     &                  nOccLeftJ(iSymJ)-nFroLeftJ(iSymJ))
            Do Lj = 1, LjOrb
               If(nFroLeftJ(iSymJ) .gt. 0) Then
                  iJ = Lj
               Else
                  iJ = iFirstS(iSymJ,jBatch) + Lj - nFro(iSymJ) - 1
               End If
               Do iB = 1,nVir(iSymB)
                  Lbj = LiPQprod(iSymB,iSymJ,jBatch)
     &                + (nFro(iSymB) + nOcc(iSymB)
     &                +  nVir(iSymB) + nDel(iSymB))
     &                * (Lj + nFroLeftJ(iSymJ) -1)
     &                + iB + nFro(iSymB) + nOcc(iSymB)
                  Do iSymI = 1 , nSym
                     iSymA = MulD2h(iSymI,iSymAI)
                     iSymAJ = MulD2h(iSymA,iSymJ)
                     iSymBI = iSymAJ
                     iSymP = iSymA
                     LiOrb = min(LnBatOrb(iSymI,iBatch)
     &                     - nFroLeftI(iSymI),
     &                       nOccleftI(iSymI)-nFroLeftI(iSymI))
                      Do Li = 1, LiOrb
                        If(nFroLeftI(iSymI) .gt. 0) Then
                           iI = Li
                        Else
                           iI  = iFirstS(iSymI,iBatch) - nFro(iSymI)
     &                      + Li - 1
                        End If
                        Lbi = LiPQprod(iSymB,iSymI,iBatch)
     &                      + (nFro(iSymB) + nOcc(iSymB)
     &                      +  nVir(iSymB) + nDel(iSymB))
     &                      * (Li + nFroLeftI(iSymI)-1)
     &                      + iB + nFro(iSymB) + nOcc(iSymB)
                        Do iA = 1, nVir(iSymA)
                           Lai = LiPQprod(iSymA,iSymI,iBatch)
     &                         + (nOcc(iSymA) + nVir(iSymA)
     &                         +  nFro(iSymA) + nDel(iSymA))
     &                         * (Li + nFroLeftI(iSymI) - 1)
     &                         + iA + nFro(iSymA) + nOcc(iSymA)
                           Laj = LiPQprod(iSymA,iSymJ,jBatch)
     &                         + (nOcc(iSymA) + nVir(iSymA)
     &                         +  nFro(iSymA) + nDel(iSymA))
     &                         * (Lj + nFroLeftJ(iSymJ) - 1)
     &                         + iA + nFro(iSymA) + nOcc(iSymA)
*
*     Construct the T_iajb (common for all terms)
*
                           Dnom1 = EOcc(iOcc(iSymI) + iI)
     &                              - EVir(iVir(iSymA) + iA)
     &                              + EOcc(iOcc(iSymJ) + iJ)
     &                              - EVir(iVir(iSymB) + iB)
                           If(iBatch.eq.jBatch) Then
                              ip_iajb = LiPQRSprod(iSymAI)
     &                                + iTri(Lai,Lbj)
                              ip_ibja = LiPQRSprod(iSymBI)
     &                                + iTri(Lbi,Laj)
                           Else
                              ip_iajb = LiPQRSprod(iSymAI)
     &                                + LnPQprod(iSymAI,iBatch)
     &                                * (Lbj-1) + Lai
                              ip_ibja = LiPQRSprod(iSymBI)
     &                                + LnPQprod(iSymBI,iBatch)
     &                                * (Laj-1) + Lbi
                              End If
                              T1 = 4.0d0*Xaibj(ip_iajb)
                              T1 = (T1 - 2.0d0*Xaibj(ip_ibja))
     &                           / Dnom1
*     Here we calculate the MP2-energy
                              EMP2_Dens = EMP2_Dens
     &                                  + 0.5d0*T1*Xaibj(ip_iajb)
                           Do iP = 1, nOrb(iSymP) + nDel(iSymP)
                              Lpi = LiPQprod(iSymP,iSymI,iBatch)
     &                            + (nOcc(iSymP) + nVir(iSymP)
     &                            +  nFro(iSymP) + nDel(iSymP))
     &                            * (Li + nFroLeftI(iSymI) - 1) + iP
                              Lpj =  LiPQprod(iSymP,iSymJ,jBatch)
     &                            + (nOcc(iSymP) + nVir(iSymP)
     &                            +  nFro(iSymP) + nDel(iSymP))
     &                            * (Lj + nFroLeftJ(iSymJ) - 1) + iP
                              If(iBatch.eq.jBatch) Then
                                 ip_ipjb = LiPQRSprod(iSymBJ)
     &                                   + iTri(Lpi,Lbj)
*                                ip_jpib = 0
                              Else
                                 ip_ipjb = LiPQRSprod(iSymBJ)
     &                                   + LnPQprod(iSymBJ,iBatch)
     &                                   * (Lbj-1) + Lpi
*                                ip_jpib = LiPQRSprod(iSymBI)
*    &                                   + LnPQprod(iSymBI,iBatch)
*    &                                   * (Lpj-1) + Lbi
                              End If
                              X2 = Xaibj(ip_ipjb)
*******************************************************************
*        Calculate P_Ca frozen virtual - active virtual
*******************************************************************
                              If(iP.gt. nFro(iSymP) + nOcc(iSymP)
     &                                + nVir(iSymP)) Then
                                 iCFroz = iP - nFro(iSymP) - nOcc(iSymP)
     &                                  - nVir(iSymP)
                                 Dnom2 = EVir(iVir(iSymA) + iA) -
     &                                   EDel(iDel(iSymP) + iCFroz)
                                 Work(iDensVirFro(iA,iCFroz,iSymP)) =
     &                                Work(iDensVirFro(iA,iCFroz,iSymP))
     &                                + T1*X2/Dnom2
                                 Work(iDensFroVir(iCFroz,iA,iSymA)) =
     &                                Work(iDensFroVir(iCFroz,iA,iSymA))
     &                                + T1*X2/Dnom2
                                 Work(iWDensVirFro(iA,iCFroz,iSymP)) =
     &                               Work(iWDensVirFro(iA,iCFroz,iSymP))
     &                               - 2.0d0 * T1*X2
***************** The 2 ON THE ROW ABOVE IS QUESTIONABLE *****'''' ***
*                 Write(6,*) 'index Dens', iDensFroVir(iCFroz,iA,iSymA) -
*     &                                    ip_Density(iSymA)
*                 Write(6,*) 'Value', Work(iDensFroVir(iCFroz,iA,iSymA))
*     &                                + 2.0d0 * T1*X2/Dnom2
*                        Write(6,*) 'xicjb', X2
*                        Write(6,*) 'EDiffbc',Dnom2
*                        Write(6,*) 'E_a', EVir(iVir(iSymA) + iA)
*                        Write(6,*) 'E:B', EDel(iDel(iSymP) + iCFroz)
*                        Write(6,*) 'Tij', T1
*                        Write(6,*) 'iB, iC',
*     &                             iCFroz+nOcc(iSymP) -nFro(iSymP), iB
*                        Write(6,*) 'iI, iJ', iI, iJ
*******************************************************************
*         Calculate P_ab active vir - active vir
*******************************************************************
                              Else If(iP .gt. nFro(iSymP) + nOcc(iSymP))
     &                             Then
                                 iC = iP - nFro(iSymP) - nOcc(iSymP)
                                 Dnom2 = + EOcc(iOcc(iSymI) + iI)
     &                                   - EVir(iVir(iSymP) + iC)
     &                                   + EOcc(iOcc(iSymJ) + iJ)
     &                                   - EVir(iVir(iSymB) + iB)
                                 Work(iDensVir(iA,iC,iSymP)) =
     &                                Work(iDensVir(iA,iC,iSymP))
     &                                + T1*X2/Dnom2
                                 Work(iWDensVir(iA,iC,iSymP)) =
     &                                Work(iWDensVir(iA,iC,iSymP))
     &                                - 2.0d0 * T1*X2
*******************************************************************
*        Calculate Lagr(3)_ai  active virtual - all occupied
*******************************************************************
                              Else
                                 iK = iP
                                Work(iMp2Lagr(iA,iK,iSymP))=
     &                                Work(iMp2Lagr(iA,iK,iSymP))
     &                                   - T1*X2
                                Work(iWDensVactOcc(iA,iK,iSymP)) =
     &                               Work(iWDensVactOcc(iA,iK,iSymP))
     &                               - 4.0d0 * T1*X2
                              End If
*
* ------------------------- Debugging comments --------------------------
*                              If(iP.gt. nFro(iSymP) + nOcc(iSymP)
*     &                             .and. (ip .le. nFro(iSymP) +
*     &                              nOcc(iSymP) + nVir(iSymP))) Then
*                                 Write(6,*) 'AIBJC',iA,iI,iB,iJ,iC
*                                 Write(6,*) 'Symm',iSymA,iSymI,iSymB,
*     &                                             iSymJ,iSymP
*                                 Write(6,*) 'Dnom1', Dnom1
*                                Write(6,*) 'Dnom2', Dnom2
*                                 Write(6,*) 'iajb', Xaibj(ip_iajb)
*                                 Write(6,*) 'ibja', Xaibj(ip_ibja)
*                                 Write(6,*) 'icjb', Xaibj(ip_ipjb)
*                                 Write(6,*) 'Bidrag T',T1*X2/dnom2
*                              Write(6,*) 'adress',iDensVir(iA,iC,iSymP)-
*     &                                    ip_Density(1)
*                              End If
*-------------------------------------------------------------------------
* ------------------------- Debugging comments --------------------------
*                              If(iP .le. nFro(iSymP) + nOcc(iSymP)) Then
*                                 Write(6,*) 'AIBJK',iA,iI,iB,iJ,iK
*                                 Write(6,*) 'Symm',iSymA,iSymI,iSymB,
*     &                                             iSymJ,iSymP
*                                 Write(6,*) 'Dnom2', Dnom2
*                                 Write(6,*) 'iajb', Xaibj(ip_iajb)
*                                 Write(6,*) 'ibja', Xaibj(ip_ibja)
*                                 Write(6,*) 'ibjk', Xaibj(ip_ipjb)
*                                 If(iBatch.ne.jBatch) Then
*                                    Write(6,*) 'ikjb', Xaibj(ip_jpib)
*                                  Write(6,*) 'Bidrag U',U1*U2
*                                 EndIf
*                                 Write(6,*) 'Bidrag T', T1*T2
*                                 Write(6,*) 'adress',
*     &                                      iMp2Lagr(iA,iK,iSymP)-
*     &                                               ip_Mp2Lagr(1)
*                              End If
* ------------------------------------------------------------------------
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do

***************************************************************************
****  Common code for Pij, Wij(I), Lagr(4) and PiJ
***************************************************************************
*
*     To avoid compiler warnings
*         iKfroz = 0
*         iK = 0

         Do iSymA = 1, nSym
            iSymI = MulD2h(iSymA,iSymAI)
            LaOrb = Min(LnBatOrb(iSymA,iBatch)-nOccLeftI(iSymA),
     &                  nOrbLeftI(iSymA)-nOccLeftI(iSymA))
             Do La = 1, LaOrb
               If(nOccLeftI(iSymA).gt.0) Then
                  iA = La
               Else
                  iA = iFirstS(iSymA,iBatch) - nFro(iSymA) - nOcc(iSymA)
     &               + La-1
               End If
               Do iI = 1, nOcc(iSymI)
                  Lai = LiPQprod(iSymI,iSymA,iBatch)
     &                + (nOcc(iSymI) + nVir(iSymI)
     &                +  nFro(iSymI) + nDel(iSymI))
     &                * (La + nOccLeftI(iSymA)-1)
     &                + iI + nFro(iSymI)
                  Do iSymJ = 1, nSym
                     iSymB = MulD2h(iSymJ,iSymBJ)
                     iSymAJ = MulD2h(iSymA,iSymJ)
                     iSymBI = iSymAJ
                     iSymP = iSymJ
                     LbOrb = Min(LnBatOrb(iSymB,jBatch)
     &                     - nOccLeftJ(iSymB),
     &                       nOrbLeftJ(iSymB)-nOccLeftJ(iSymB))
                     Do Lb = 1, LbOrb
                        If(nOccLeftJ(iSymB).gt.0) Then
                           iB = Lb
                        Else
                           iB = iFirstS(iSymB,jBatch)
     &                        - nFro(iSymB) - nOcc(iSymB) + Lb-1
                        End If
                        Lbi = LiPQprod(iSymI,iSymB,jBatch)
     &                      + (nOcc(iSymI) + nVir(iSymI)
     &                      +  nFro(iSymI) + nDel(iSymI))
     &                      * (Lb+nOccLeftJ(iSymB)-1)
     &                      + iI + nFro(iSymI)
                        Do iJ = 1, nOcc(iSymJ)
                           Lbj = LiPQprod(iSymJ,iSymB,jBatch)
     &                         + (nOcc(iSymJ) + nVir(iSymJ)
     &                         +  nFro(iSymJ) + nDel(iSymJ))
     &                         * (Lb+nOccLeftJ(iSymB)-1)
     &                         + iJ + nFro(iSymJ)
                           Laj = LiPQprod(iSymJ,iSymA,iBatch)
     &                         + (nOcc(iSymJ) + nVir(iSymJ)
     &                         +  nFro(iSymJ) + nDel(iSymJ))
     &                         * (La+nOccLeftI(iSymA)-1)
     &                         + iJ + nFro(iSymJ)
                           Dnom1 = EOcc(iOcc(iSymI) + iI)
     &                           - EVir(iVir(iSymA) + iA)
     &                           + EOcc(iOcc(iSymJ) + iJ)
     &                           - EVir(iVir(iSymB) + iB)
                           If(iBatch.eq.jBatch) Then
                              ip_iajb = LiPQRSprod(iSymAI)
     &                                + iTri(Lai,Lbj)
                              ip_ibja = LiPQRSprod(iSymBI)
     &                                + iTri(Lbi,Laj)
                           Else
                              ip_iajb = LiPQRSprod(iSymAI)
     &                                + LnPQprod(iSymAI,iBatch)
     &                                * (Lbj-1) + Lai
                              ip_ibja = LiPQRSprod(iSymBI)
     &                                + LnPQprod(iSymBI,iBatch)
     &                                * (Lbi-1) + Laj
                           End If
                           T1 = 4.0D0*Xaibj(ip_iajb)
                           T1 = (T1-2.0D0*Xaibj(ip_ibja))
     &                            / Dnom1
                           Do iP = 1, nOrb(iSymP) + nDel(iSymP)
*
                              Lbp = LiPQprod(iSymP,iSymB,jBatch)
     &                            + (nOcc(iSymP) + nVir(iSymP)
     &                            +  nFro(iSymP) + nDel(iSymP))
     &                            * (Lb+nOccLeftJ(iSymB)-1)
     &                            + iP
                              If(iBatch.eq.jBatch) Then
                                 ip_iapb = LiPQRSprod(iSymAI)
     &                                   + iTri(Lai,Lbp)
                              Else
                                 ip_iapb = LiPQRSprod(iSymAI)
     &                                   + LnPQprod(iSymAI,iBatch)
     &                                   * (Lbp-1) + Lai
                              End If
*
                              X2 = Xaibj(ip_iapb)
*******************************************************************
*        Calculate P_Pj frozen occupied - active occupied
*******************************************************************
                              If(ip .le. nFro(iSymP)) Then
                                 iKFroz = iP
                                 Dnom2 = EOcc(iOcc(iSymI) + iI)
     &                                 - EFro(iFro(iSymP) + iKFroz)
                                 Work(iDensFroOcc(iKFroz, iJ,iSymJ)) =
     &                             Work(iDensFroOcc(iKFroz,iJ,iSymJ))
     &                             + T1*X2/Dnom2
                                 Work(iDensOccFro(iJ,iKFroz,iSymJ)) =
     &                             Work(iDensOccFro(iJ,iKFroz,iSymJ))
     &                             + T1*X2/Dnom2
                                 Work(iWDensOccFro(iJ,iKFroz,iSymP)) =
     &                             Work(iWDensOccFro(iJ,iKFroz,iSymP))
     &                             - 2.0d0 * T1*X2
*******************************************************************
*         Calculate P_ij active occ - active occ
*******************************************************************
                              Else If(iP .le. nFro(iSymP)
     &                                      + nOcc(iSymP)) Then
                                 iK = iP - nFro(iSymP)
                                 Dnom2 = EOcc(iOcc(iSymI) + iI)
     &                                  - EVir(iVir(iSymA) + iA)
     &                                  + EOcc(iOcc(iSymP) + iK)
     &                                  - EVir(iVir(iSymB) + iB)
                                 Work(iDensOcc(iJ,iK,iSymP)) =
     &                                Work(iDensOcc(iJ,iK,iSymP))
     &                                - T1*X2/Dnom2
                                 Work(iWDensOcc(iJ,iK,iSymP)) =
     &                             Work(iWDensOcc(iJ,iK,iSymP))
     &                             - 2.0d0 * T1*X2
*******************************************************************
*        Calculate Lagr(4)_ai all virtual - active occupied
*******************************************************************
                              Else
                                 iC = iP - nFro(iSymP) - nOcc(iSymP)
                                 Work(iMp2Lagr(iC,iJ+nFro(iSymJ),
     &                                iSymJ)) =
     &                           Work(iMp2Lagr(iC,iJ+nFro(iSymJ),iSymJ))
     &                                + T1*X2
                              End If

* ------------------------- Debugging comments --------------------------
*                           If((iP .gt. nFro(iSymP)).and.
*     &                        (iP .le. nFro(iSymP)+nOcc(iSymP))) Then
*                              Write(6,*) 'AIBJK',iA,iI,iB,iJ,iK
*                              Write(6,*) 'Symm',iSymA,iSymI,iSymB,
*     &                                          iSymJ,iSymP
*                              Write(6,*) 'Dnom1', Dnom1
*                              Write(6,*) 'Dnom2', Dnom2
*                              Write(6,*) 'iajb', Xaibj(ip_iajb)
*                              Write(6,*) 'ibja', Xaibj(ip_ibja)
*                              Write(6,*) 'iakb', Xaibj(ip_iapb)
*                              Write(6,*) 'Bidrag', T1*T2
*                              Write(6,*) 'adress',iDensOcc(iJ,iK,iSymP)-
*     &                                            ip_Density(1)
*                           End If
* ----------------------------------------------------------------------
* ------------------------- Debugging comments -------------------------
*                              Write(6,*) 'AIBJC',iA,iI,iB,iJ,iC
*                              Write(6,*) 'Symm',iSymA,iSymI,iSymB,
*     &                                          iSymJ,iSymC
*                              Write(6,*) 'Dnom1', Dnom1
*                              Write(6,*) 'aibj', Xaibj(ip_iajb)
*                              Write(6,*) 'biaj', Xaibj(ip_ibja)
*                              Write(6,*) 'iacb', Xaibj(ip_iacb)
*                              If(iBatch.ne.jBatch) Then
*                                 Write(6,*) 'iacb', Xaibj(ip_ibca)
*                                 Write(6,*) 'Bidrag U', U1*U2
*                              End If
*                              Write(6,*) 'Bidrag T', T1*T2
*                              Write(6,*) 'adress',iMp2Lagr(iC,iJ,iSymJ)-
*     &                                   ip_Mp2Lagr(1)
*------------------------------------------------------------------------
*
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      Return
      End
