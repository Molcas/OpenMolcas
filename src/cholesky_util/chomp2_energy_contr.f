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
* Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
*               Francesco Aquilante                                    *
************************************************************************
      SubRoutine ChoMP2_Energy_Contr(EMP2,EOcc,EVir,Xaibj,LnT2am,LiT2am,
     &                               iBatch,jBatch)
C
C     Thomas Bondo Pedersen, Dec. 2004 / Feb. 2005.
C
C     Purpose: compute (MINUS the) energy contribution from a
C              batch of (ai|bj) integrals.
C
C     Modified by F. Aquilante to compute separately the Opposite-Spin
C                              contribution to the MP2 energy
#include "implicit.fh"
      Real*8  EOcc(*), EVir(*), Xaibj(LnT2am)
      Integer LiT2am(8)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Integer a, b, aibj, biaj, ab, ba, abij, baij

      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      iFirstS(i,j)=iWork(ip_FirstS-1+nSym*(j-1)+i)
      LnOcc(i,j)=iWork(ip_LnOcc-1+nSym*(j-1)+i)
      LiT1am(i,j,k)=iWork(ip_LiT1am-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      LnT1am(i,j)=iWork(ip_LnT1am-1+nSym*(j-1)+i)
      LiMatij(i,j,k)=iWork(ip_LiMatij-1+nSym*nSym*(k-1)+nSym*(j-1)+i)

      If (DoT1amp) Then
         Call ChoMP2_Energy_Contr_T1(EMP2,EOcc,EVir,Xaibj,LnT2am,LiT2am,
     &                               iBatch,jBatch)
         Return
      EndIf

      Call Cho_GAdGOp(Xaibj,LnT2am,'+')

      If (iBatch .eq. jBatch) Then

         If (ChoAlg .eq. 2) Then ! M(ab,ij)=(ai|bj) with i<=j.
            Do iSymj = 1,nSym
               iSymi = iSymj
               Do Lj = 1,LnOcc(iSymj,jBatch)
                  j = iFirstS(iSymj,jBatch) + Lj - 1
                  Do Li = 1,LnOcc(iSymi,iBatch)
                     i = iFirstS(iSymi,iBatch) + Li - 1
                     ij = LiMatij(iSymi,iSymj,iBatch) + iTri(Li,Lj)
                     Do iSymb = 1,nSym
                        iSyma = iSymb
                        Do b = 1,nVir(iSymb)
                           Do a = 1,nVir(iSyma)
                              ab = iMatab(iSyma,iSymb)
     &                           + nVir(iSyma)*(b-1) + a
                              ba = iMatab(iSymb,iSyma)
     &                           + nVir(iSymb)*(a-1) + b
                              abij = LiT2am(1)
     &                             + nMatab(1)*(ij-1)
     &                             + ab
                              baij = LiT2am(1)
     &                             + nMatab(1)*(ij-1)
     &                             + ba
                              Dnom = EVir(iVir(iSyma)+a)
     &                             - EOcc(iOcc(iSymi)+i)
     &                             + EVir(iVir(iSymb)+b)
     &                             - EOcc(iOcc(iSymj)+j)
                              Taibj = Xaibj(abij)/Dnom
                              Waibj = 2.0D0*Xaibj(abij)
                              EOSMP2= EOSMP2 + Taibj*Waibj
                              Waibj = Waibj - Xaibj(baij)
                              Eaibj = Taibj*Waibj
                              EMP2  = EMP2 + Eaibj
                              WREF  = WREF + Eaibj/Dnom
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
            Do iSymij = 2,nSym
               iSymab = iSymij
               Do iSym2 = 1,nSym
                  iSym1 = MulD2h(iSym2,iSymij)
                  iSymj = max(iSym1,iSym2)
                  iSymi = min(iSym1,iSym2)
                  Do Lj = 1,LnOcc(iSymj,jBatch)
                     j = iFirstS(iSymj,jBatch) + Lj - 1
                     Do Li = 1,LnOcc(iSymi,iBatch)
                        i = iFirstS(iSymi,iBatch) + Li - 1
                        ij = LiMatij(iSymi,iSymj,iBatch)
     &                     + LnOcc(iSymi,iBatch)*(Lj-1) + Li
                        Do iSymb = 1,nSym
                           iSyma = MulD2h(iSymb,iSymab)
                           Do b = 1,nVir(iSymb)
                              Do a = 1,nVir(iSyma)
                                 ab = iMatab(iSyma,iSymb)
     &                              + nVir(iSyma)*(b-1) + a
                                 ba = iMatab(iSymb,iSyma)
     &                              + nVir(iSymb)*(a-1) + b
                                 abij = LiT2am(iSymij)
     &                                + nMatab(iSymab)*(ij-1)
     &                                + ab
                                 baij = LiT2am(iSymij)
     &                                + nMatab(iSymab)*(ij-1)
     &                                + ba
                                 Dnom = EVir(iVir(iSyma)+a)
     &                                - EOcc(iOcc(iSymi)+i)
     &                                + EVir(iVir(iSymb)+b)
     &                                - EOcc(iOcc(iSymj)+j)
                                 Taibj = Xaibj(abij)/Dnom
                                 Waibj = 2.0D0*Xaibj(abij)
                                 EOSMP2= EOSMP2 + Taibj*Waibj
                                 Waibj = Waibj - Xaibj(baij)
                                 Eaibj = Taibj*Waibj
                                 EMP2  = EMP2 + Eaibj
                                 WREF  = WREF + Eaibj/Dnom
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         Else ! triangular storage, (ai|bj) with ai<=bj.
            Do iSymbj = 1,nSym
               iSymai = iSymbj
               Do iSymj = 1,nSym
                  iSymb = MulD2h(iSymj,iSymbj)
                  Do Lj = 1,LnOcc(iSymj,jBatch)
                     j = iFirstS(iSymj,jBatch) + Lj - 1
                     Do b = 1,nVir(iSymb)
                        Lbj = LiT1am(iSymb,iSymj,jBatch)
     &                      + nVir(iSymb)*(Lj - 1) + b
                        Do iSymi = 1,nSym
                           iSyma  = MulD2h(iSymi,iSymai)
                           iSymaj = MulD2h(iSyma,iSymj)
                           iSymbi = iSymaj
                           Do Li = 1,LnOcc(iSymi,iBatch)
                              i   = iFirstS(iSymi,iBatch) + Li - 1
                              Lbi = LiT1am(iSymb,iSymi,iBatch)
     &                            + nVir(iSymb)*(Li - 1) + b
                              Do a = 1,nVir(iSyma)
                                 Lai = LiT1am(iSyma,iSymi,iBatch)
     &                               + nVir(iSyma)*(Li - 1) + a
                                 Laj = LiT1am(iSyma,iSymj,jBatch)
     &                               + nVir(iSyma)*(Lj - 1) + a
                                 aibj = LiT2am(iSymai) + iTri(Lai,Lbj)
                                 biaj = LiT2am(iSymbi) + iTri(Lbi,Laj)
                                 Dnom = EVir(iVir(iSyma)+a)
     &                                - EOcc(iOcc(iSymi)+i)
     &                                + EVir(iVir(iSymb)+b)
     &                                - EOcc(iOcc(iSymj)+j)
                                 Taibj = Xaibj(aibj)/Dnom
                                 Waibj = 2.0D0*Xaibj(aibj)
                                 EOSMP2= EOSMP2 + Taibj*Waibj
                                 Waibj = Waibj - Xaibj(biaj)
                                 Eaibj = Taibj*Waibj
                                 EMP2  = EMP2 + Eaibj
                                 WREF  = WREF + Eaibj/Dnom
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End If

      Else ! rectangular storage (ai|bj) with ai<bj.

         EMP2_sav = EMP2
         EMP2     = 0.0D0
         WREF_sav = WREF
         WREF     = 0.0D0
         EOSMP2_sav = EOSMP2
         EOSMP2     = 0.0D0
         Do iSymbj = 1,nSym
            iSymai = iSymbj
            Do iSymj = 1,nSym
               iSymb = MulD2h(iSymj,iSymbj)
               Do Lj = 1,LnOcc(iSymj,jBatch)
                  j = iFirstS(iSymj,jBatch) + Lj - 1
                  Do b = 1,nVir(iSymb)
                     Lbj = LiT1am(iSymb,iSymj,jBatch)
     &                   + nVir(iSymb)*(Lj - 1) + b
                     Do iSymi = 1,nSym
                        iSyma  = MulD2h(iSymi,iSymai)
                        iSymaj = MulD2h(iSyma,iSymj)
                        iSymbi = iSymaj
                        Do Li = 1,LnOcc(iSymi,iBatch)
                           i   = iFirstS(iSymi,iBatch) + Li - 1
                           Lbi = LiT1am(iSymb,iSymi,iBatch)
     &                         + nVir(iSymb)*(Li - 1) + b
                           Do a = 1,nVir(iSyma)
                              Lai = LiT1am(iSyma,iSymi,iBatch)
     &                            + nVir(iSyma)*(Li - 1) + a
                              Laj = LiT1am(iSyma,iSymj,jBatch)
     &                            + nVir(iSyma)*(Lj - 1) + a
                              aibj = LiT2am(iSymai)
     &                             + LnT1am(iSymai,iBatch)*(Lbj-1) + Lai
                              biaj = LiT2am(iSymbi)
     &                             + LnT1am(iSymbi,iBatch)*(Laj-1) + Lbi
                              Dnom = EVir(iVir(iSyma)+a)
     &                             - EOcc(iOcc(iSymi)+i)
     &                             + EVir(iVir(iSymb)+b)
     &                             - EOcc(iOcc(iSymj)+j)
                              Taibj = Xaibj(aibj)/Dnom
                              Waibj = 2.0D0*Xaibj(aibj)
                              EOSMP2= EOSMP2 + Taibj*Waibj
                              Waibj = Waibj - Xaibj(biaj)
                              Eaibj = Taibj*Waibj
                              EMP2  = EMP2 + Eaibj
                              WREF  = WREF + Eaibj/Dnom
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
         EMP2 = EMP2_sav + 2.0D0*EMP2
         EOSMP2 = EOSMP2_sav + 2.0D0*EOSMP2
         WREF = WREF_sav + 2.0D0*WREF

      End If

      EOSMP2 = 0.5D0*EOSMP2

      End
