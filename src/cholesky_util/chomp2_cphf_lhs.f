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

      SubRoutine ChoMP2_cphf_lhs(EOcc,EVir,EFro,EDel,Xaibj,
     &                          LnPQRSprod,LiPQRSprod,iBatch,jBatch,
     &                          ip_Ap,ip_P,nVec,iVecOff,
     &                          nOccLeftI,nOccLeftJ,
     &                          nOrbLeftI, nOrbLeftJ)
*     This will calculate a trial lefthandside of the cphf equations.
*
*
#include "implicit.fh"
      Real*8 EOcc(*), Evir(*), EFro(*), EDel(*), Xaibj(LnPQRSprod)
      Integer LiPQRSprod(8)
      Integer iVecOff(8)
      Integer nOccLeftI(8),nOccLeftJ(8)
      Integer nOrbLeftI(8),nOrbLeftJ(8)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
*
      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      iFirstS(i,j)=iWork(ip_FirstS-1+nSym*(j-1)+i)
      iAp(i,j,k)= ip_Ap + iVecOff(k) + j-1 + (nOcc(k)+nFro(k))*(i-1)
      iPvec(i,j,k)=ip_P + iVecOff(k) + j-1 + (nOcc(k)+nFro(k))*(i-1)
      LnBatOrb(i,j)=iWork(ip_LnBatOrb-1+nSym*(j-1)+i)
      LiPQprod(i,j,k) = iWork(ip_LiPQprod-1
     &                      + nSym*nSym*(k-1)+nSym*(j-1)+i)
      LnPQprod(i,j)=iWork(ip_LnPQprod-1+nSym*(j-1)+i)
*
      Do iSymJ = 1, nSym
         iSymB = iSymJ
         iSymBJ = MulD2h(iSymJ,iSymB)
         iSymAI = iSymBJ
         LjOrb = Min(LnBatOrb(iSymJ,jBatch),nOccLeftJ(iSymJ))
         Do Lj = 1,LjOrb
            iJ = iFirstS(iSymJ,jBatch) + Lj - 1
            Do iB = 1, nVir(iSymB) + nDel(iSymB)
               Lbj = LiPQprod(iSymB,iSymJ,jBatch)
     &             + (nFro(iSymB) + nOcc(iSymB)
     &             +  nVir(iSymB) + nDel(iSymB))
     &             * (Lj-1)
     &             + iB + nFro(iSymB) +  nOcc(iSymB)
               Do iSymI = 1, nSym
                  iSymA = iSymI
                  iSymAJ = MulD2h(iSymA,iSymJ)
                  iSymBI = iSymAJ
                  LiOrb = min(LnBatOrb(iSymI,iBatch),
     &                        nOccLeftI(iSymI))
                  Do Li = 1, LiOrb
                     iI  = iFirstS(iSymI,iBatch) + Li - 1
                     Lbi = LiPQprod(iSymB,iSymI,iBatch)
     &                   + (nFro(iSymB) + nOcc(iSymB)
     &                   +  nVir(iSymB) + nDel(iSymB))
     &                   * (Li-1)
     &                   + iB + nFro(iSymB) + nOcc(iSymB)
                     Do iA = 1, nVir(iSymA) + nDel(iSymA)
                        Lai = LiPQprod(iSymA,iSymI,iBatch)
     &                      + (nFro(iSymA) + nOcc(iSymA)
     &                      +  nVir(iSymA) + nDel(iSymA))
     &                      * (Li - 1)
     &                      + iA + nFro(iSymA) + nOcc(iSymA)
                        Laj = LiPQprod(iSymA,iSymJ,jBatch)
     &                      + (nFro(iSymA) + nOcc(iSymA)
     &                      +  nVir(iSymA) + nDel(iSymA))
     &                      * (Lj - 1)
     &                      + iA + nFro(iSymA) + nOcc(iSymA)
*
                        If(iBatch.eq.jBatch) Then
                           ip_aibj = LiPQRSprod(iSymAI)
     &                             + iTri(Lai,Lbj)
                           ip_biaj = LiPQRSprod(iSymBI)
     &                             + iTri(Lbi,Laj)
                        Else
                           ip_aibj = LiPQRSprod(iSymAI)
     &                             + LnPQprod(iSymAI,iBatch)
     &                             * (Lbj-1) + Lai
                           ip_biaj = LiPQRSprod(iSymBI)
     &                             + LnPQprod(iSymBI,iBatch)
     &                             * (Laj-1) + Lbi
                        End If
                        Diag1 = 4.0d0 * Xaibj(ip_aibj)
     &                        -         Xaibj(ip_biaj)
*
                        If((iI.eq.iJ).and.(iB.eq.iA) .and.
     &                     (iSymI.eq.iSymJ))         Then
                           If(iA .le. nVir(iSymA)) Then
                              E_Vir = EVir(iVir(iSymA) + iA)
                           Else
                              E_Vir = EDel(iDel(iSymA) + iA-nVir(iSymA))
                           End If
                           If(iI .gt. nFro(iSymI)) Then
                              E_Occ = EOcc(iOcc(iSymI) + iI-nFro(iSymI))
                           Else
                              E_Occ = EFro(iFro(iSymI) + iI)
                           End If
                           Ediff = E_Vir - E_Occ
                           Diag1 = Diag1 + Ediff
                        End If
*
                        Work(iAp(iB,iJ,iSymJ)) =
     &                       Work(iAp(iB,iJ,iSymJ)) +
     &                       Diag1*Work(iPvec(iA,iI,iSymI))
*-------------- Debug comments ---------------------------------------------
*                        Write(6,*) 'Symm A B', iSymA, iSymB
*                        Write(6,*) 'IAJB', iI,iA,iJ,iB
*                        Write(6,*) 'aibj', Xaibj(ip_aibj)
*                        Write(6,*) 'biaj', Xaibj(ip_biaj)
*                        Write(6,*) 'Ediff', Ediff
*                        Write(6,*) 'contr exch + Energy',
*     &                              Diag1*Work(iPvec(iA,iI,iSymI))
*                        Write(6,*) 'P',Work(iPVec(iA,iI,iSymI))
*                        Write(6,*) 'adress AP',
*     &                              iAp(iB,iJ,iSymJ)-ip_AP
*-----------------------------------------------------------------------------
                     End Do
                  End Do
               End Do
            End Do
*****************************************************************************
****     Coulomb part
*****************************************************************************
            Do iSymI = 1, nSym
               iSymA = iSymI
               iSymAJ = MulD2h(iSymA,iSymJ)
               iSymIJ = MulD2h(iSymI,iSymJ)
               Do iI = 1, nOcc(iSymI) + nFro(iSymI)
                  Lij = LiPQprod(iSymI,iSymJ,jBatch)
     &                + (nFro(iSymI) + nOcc(iSymI)
     &                +  nVir(iSymI) + nDel(iSymI))
     &                * (Lj-1) + iI
                  LbOrb = Min(LnBatOrb(iSymB,iBatch)
     &                     -  nOccLeftI(iSymB),
     &                        nOrbLeftI(iSymB)-nOccLeftI(iSymB))
                  Do Lb = 1, LbOrb
                     If(nOccLeftI(iSymB).gt.0) Then
                        iB = Lb
                     Else
                        iB = iFirstS(iSymB,iBatch)
     &                     - nOcc(iSymB) - nFro(iSymB) + Lb-1
                     End If
                     Do iA = 1, nVir(iSymA) + nDel(iSymA)
                        Lab = LiPQprod(iSymA,iSymB,iBatch)
     &                      + (nFro(iSymA) + nOcc(iSymA)
     &                      +  nVir(iSymA) + nDel(iSymA))
     &                      * (Lb + nOccLeftI(iSymB) - 1)
     &                      + iA + nFro(iSymA) + nOcc(iSymA)
*
                        If(iBatch.eq.jBatch) Then
                           ip_jiba = LiPQRSprod(iSymIJ)
     &                             + iTri(Lij,Lab)
                        Else
                           ip_jiba = LiPQRSprod(iSymIJ)
     &                             + LnPQprod(iSymIJ,iBatch)
     &                             * (Lij - 1) + Lab
                        End If
*
                        Diag2 = Xaibj(ip_jiba)
                        Work(iAp(iB,iJ,iSymJ)) =
     &                       Work(iAp(iB,iJ,iSymJ)) -
     &                       Diag2*Work(iPvec(iA,iI,iSymI))
*------------ Debug comments -------------------------------------------
*                        Write(6,*) 'Symm A B', iSymA, iSymB
*                        Write(6,*) 'IAJB', iI,iA,iJ,iB
*                        Write(6,*) 'jiba', Diag2
*                        Write(6,*) 'adress AP',
*     &                              iAp(iB,iJ,iSymJ)-ip_AP
*                        Write(6,*) 'Contr coul',
*     &                              Diag2*Work(iPvec(iA,iI,iSymI))
*------------------------------------------------------------------------
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nVec)
         Call Unused_integer_array(nOrbLeftJ)
      End If
      End
