!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!***********************************************************************
      SubRoutine ChoMP2_Energy_GetInd(LnT2am,LiT2am,iBatch,jBatch)
!
!     Thomas Bondo Pedersen, Dec. 2004 / Feb. 2005.
!
!     Purpose: setup (ai|bj) index arrays for batch i,j.
!              For iBatch=jBatch and ChoAlg=2, (ai|bj) is stored as
!              the matrix M(ab,ij) with i<=j.
!
      use ChoMP2, only: LnT1am, LnMatij
      Implicit None
      Integer LnT2am, iBatch, jBatch
      Integer LiT2am(8)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"

      Integer iSym

      Character*14 String
      Character*20 SecNam
      Parameter (SecNam = 'ChoMP2_Energy_GetInd')

      If (iBatch .eq. jBatch) Then
         LnT2am = 0
         If (ChoAlg .eq. 1) Then
            Do iSym = 1,nSym
               LiT2am(iSym) = LnT2am
               LnT2am = LnT2am
     &                + LnT1am(iSym,iBatch)*(LnT1am(iSym,iBatch)+1)/2
            End Do
         Else If (ChoAlg .eq. 2) Then
            Do iSym = 1,nSym
               LiT2am(iSym) = LnT2am
               LnT2am = LnT2am
     &                + nMatab(iSym)*LnMatij(iSym,iBatch)
            End Do
         Else
            Write(String,'(A8,I6)') 'ChoAlg =',ChoAlg
            Call ChoMP2_Quit(SecNam,'ChoAlg out-of-bounds error!',
     &                       String)
         End If
      Else
         LnT2am = 0
         Do iSym = 1,nSym
            LiT2am(iSym) = LnT2am
            LnT2am = LnT2am + LnT1am(iSym,iBatch)*LnT1am(iSym,jBatch)
         End Do
      End If

      End
