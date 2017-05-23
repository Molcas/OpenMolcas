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
      SubRoutine ChoMP2_Energy_GetPQInd(LnPQRSprod, LiPQRSprod,
     &                                  iBatch,jBatch)
C
C     Jonas Bostrom, june 2008
C
C     Purpose: setup (pq|rs) index arrays for calculating mp2_densities.
C
      Implicit None
      Integer iBatch, jBatch
      Integer LnPQRSprod,LiPQRSprod(8)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Integer iSym
      Integer i, j
      Integer LnPQprod

      Character*14 String
      Character*22 SecNam
      Parameter (SecNam = 'ChoMP2_Energy_GetPQInd')

      LnPQprod(i,j)=iWork(ip_LnPQprod-1+nSym*(j-1)+i)

      If (iBatch .eq. jBatch) Then
         LnPQRSprod = 0
         If (ChoAlg .eq. 1) Then
            Do iSym = 1,nSym
               LiPQRSprod(iSym) = LnPQRSprod
               LnPQRSprod = LnPQRSprod
     &                    + LnPQprod(iSym,iBatch)
     &                    * (LnPQprod(iSym,iBatch)+1)/2
            End Do
         Else
            Write(String,'(A8,I6)') 'ChoAlg =',ChoAlg
            Call qEnter(SecNam)
            Call ChoMP2_Quit(SecNam,'ChoAlg out-of-bounds error!',
     &                       String)
         End If
      Else
         LnPQRSprod = 0
         Do iSym = 1,nSym
            LiPQRSprod(iSym) = LnPQRSprod
            LnPQRSprod = LnPQRSprod + LnPQprod(iSym,iBatch)*
     &                                LnPQprod(iSym,jBatch)
         End Do
      End If

      End
