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
      SubRoutine Cho_P_SetGL()
C
C     Purpose: set global and local index arrays and diagonal.
C              If a sequencial run:
C                 Diag   => Diag_Hidden
C                 Diag_G => Null()
C              If a parallel run:
C                 Diag   => Diag_G_Hidden
C                 Diag_G => Diag_Hidden
C
C              Diag_Hidden is allocated in the calling routine, while
C              Diag_G_Hidden is allocated here if needed.
C
      use ChoSwp, only: nnBstRSh, nnBstRSh_G, nnBstRsh_L_Hidden
      use ChoSwp, only: iiBstRSh, iiBstRSh_G, iiBstRsh_L_Hidden
      use ChoSwp, only: IndRSh, IndRSh_G, IndRsh_G_Hidden
      use ChoSwp, only: InfRed, InfRed_G, InfRed_G_Hidden
      use ChoSwp, only: InfVec, InfVec_G, InfVec_G_Hidden
      use ChoSwp, only: IndRed, IndRed_G, IndRed_G_Hidden
      use ChoSwp, only: Diag, Diag_G, Diag_Hidden, Diag_G_Hidden
      use ChoArr, only: iL2G
      Implicit None
#include "cholesky.fh"
#include "choptr2.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Character*11 SecNam
      Parameter (SecNam = 'Cho_P_SetGL')

      Integer N, iSP, iSym, iShlAB, i1, i2, irc

      Integer i
      Integer mySP

      mySP(i)=iWork(ip_mySP-1+i)

C     If not parallel, return.
C     ------------------------

      If (.not.Cho_Real_Par) Then
         Diag => Diag_Hidden
         Return
      End If

C     Set global data (choglob.fh).
C     ------------------------------

      Diag_G => Diag_Hidden

      nnShl_G = nnShl
      mmBstRT_G = mmBstRT

      N = 8*3
      Call iCopy(N,iiBstR,1,iiBstR_G,1)
      Call iCopy(N,nnBstR,1,nnBstR_G,1)
      Call iCopy(3,nnBstRT,1,nnBstRT_G,1)

      InfRed_G => InfRed

      InfVec_G => InfVec

      iiBstRSh_G => iiBstRSh

      nnBstRSh_G => nnBstRSh

      IndRed_G => IndRed

      IndRSh_G => IndRSh

C     Reallocate and reset local data.
C     --------------------------------

      Call mma_allocate(InfRed_G_Hidden,SIZE(InfRed),
     &                  Label='InfRed_G_Hidden')
      InfRed => InfRed_G_Hidden

      Call mma_allocate(InfVec_G_Hidden,SIZE(InfVec,1),SIZE(InfVec,2),
     &                  SIZE(InfVec,3),Label='InfVec_G_Hidden')
      InfVec => InfVec_G_Hidden

      nnShl = n_mySP
      Call mma_allocate(iiBstRsh_L_Hidden,nSym,n_mySP,3,
     &                  Label='iiBstRSh_L_Hidden')
      iiBstRSh => iiBstRSh_L_Hidden
      Call mma_allocate(nnBstRsh_L_Hidden,nSym,n_mySP,3,
     &                  Label='nnBstRSh_L_Hidden')
      nnBstRSh => nnBstRSh_L_Hidden

      Do iSP = 1,nnShl
         iShlAB = mySP(iSP)
         Do iSym = 1,nSym
            nnBstRSh(iSym,iSP,1) = nnBstRSh_G(iSym,iShlAB,1)
         End Do
      End Do
      Call Cho_SetRedInd(1)
      mmBstRT = nnBstRT(1)

      Call mma_allocate(IndRed_G_Hidden,mmBstRT,3,
     &                  Label='IndRed_G_Hidden')
      IndRed => IndRed_G_Hidden
      Call mma_allocate(IndRSh_G_Hidden,mmBstRT,
     &                  Label='IndRSh_G_Hidden')
      IndRSh => IndRSh_G_Hidden
      call mma_allocate(iL2G,mmBstRT,Label='iL2G')

      N = 0
      Do iSym = 1,nSym
         Do iSP = 1,nnShl
            iShlAB = mySP(iSP)
            i1 = iiBstR_G(iSym,1) + iiBstRSh_G(iSym,iShlAB,1) + 1
            i2 = i1 + nnBstRSh_G(iSym,iShlAB,1) - 1
            Do i = i1,i2
               N = N + 1
               IndRed(N,1) = IndRed_G(i,1)
               IndRSh(N) = IndRSh_G(i)
               iL2G(N) = i
            End Do
         End Do
      End Do
      Call Cho_X_RSCopy(irc,1,2)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': [1] Cho_X_RSCopy returned ',irc
         Call Cho_Quit('Error in '//SecNam,104)
      End If
      Call Cho_X_RSCopy(irc,2,3)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': [2] Cho_X_RSCopy returned ',irc
         Call Cho_Quit('Error in '//SecNam,104)
      End If

C     Allocate and set local diagonal.
C     --------------------------------

      Call mma_allocate(Diag_G_Hidden,mmBstRT,Label='Diag_G_Hidden')
      Diag => Diag_G_Hidden
      Do i = 1,mmBstRT
         Diag(i) = Diag_G(iL2G(i))
      End Do

      End
