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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Srt(Vec,Srt,nVec,iSym,iBatch)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: copy out subblock of vectors.
C
#include "implicit.fh"
      Real*8  Vec(*), Srt(*)
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
#include "chomp2_cfg.fh"

      MulD2h(i,j)=iEor(i-1,j-1)+1
      LnOcc(i,j)=iWork(ip_LnOcc-1+nSym*(j-1)+i)
      LnBatOrb(i,j)=iWork(ip_LnBatOrb-1+nSym*(j-1)+i)
      LnT1am(i,j)=iWork(ip_LnT1am-1+nSym*(j-1)+i)
      LnPQprod(i,j)=iWork(ip_LnPQprod-1+nSym*(j-1)+i)
      LiT1am(i,j,k)=iWork(ip_LiT1am-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      LiPQprod(i,j,k)=iWork(ip_LiPQprod-1+nSym*nSym*(k-1)+nSym*(j-1)+i)
      iFirstS(i,j)=iWork(ip_FirstS-1+nSym*(j-1)+i)

      If(.not.DoDens) Then
         Do iVec = 1,nVec

            kOff0 = nT1am(iSym)*(iVec-1) + 1
            kOff1 = LnT1am(iSym,iBatch)*(iVec-1) + 1

            Do iSymi = 1,nSym

               iSyma = MulD2h(iSymi,iSym)
               If (LnOcc(iSymi,iBatch).gt.0 .and. nVir(iSyma).gt.0) Then
                  lCp = nVir(iSyma)*LnOcc(iSymi,iBatch)
                  kOff2 = kOff0 + iT1am(iSyma,iSymi)
     &                  + nVir(iSyma)*(iFirstS(iSymi,iBatch)-1)
                  kOff3 = kOff1 + LiT1am(iSyma,iSymi,iBatch)
                  Call dCopy_(lCp,Vec(kOff2),1,Srt(kOff3),1)
               End If

            End Do

         End Do
      Else
*       Special sorting for Mp2-density calculations where all integrals
*       are used (for pure energy calculations only integrals of type
*       (occ,vir|occ,vir) are needed).
         Do iVec = 1,nVec
*
            kOff0 = nPQ_prod(iSym)*(iVec-1) + 1
            kOff1 = LnPQprod(iSym,iBatch)*(iVec-1) + 1
*
            Do iSymI = 1,nSym
*
               iSymA = MulD2h(iSymI,iSym)
               If (LnBatOrb(iSymI,iBatch).gt.0 .and.
     &             (nFro(iSymA) + nOcc(iSymA)
     &            + nVir(iSymA) + nDel(iSymA)).gt.0) Then
                  lCp = (nFro(iSymA) + nOcc(iSymA)
     &                +  nVir(iSymA) + nDel(iSymA))
     &                *  LnBatOrb(iSymI,iBatch)
                  kOff2 = kOff0 + iPQ_prod(iSymA,iSymI)
     &                  + (nFro(iSymA) + nOcc(iSymA)
     &                  +  nVir(iSymA) + nDel(iSymA))
     &                  * (iFirstS(iSymi,iBatch)-1)
                  kOff3 = kOff1 + LiPQprod(iSyma,iSymi,iBatch)
                  Call dCopy_(lCp,Vec(kOff2),1,Srt(kOff3),1)
               End If
            End Do
         End Do
      End If


      End
