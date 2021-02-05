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
* Copyright (C) 2007, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_P_Distrib_SP(iOpt,mySP,N_mySP)
C
C     Thomas Bondo Pedersen, June 2007.
C
C     Determine distribution of ShellPairs.
C     iOpt=1: each node has same number of ShellPairs.
C     iOpt=2: each node has same dimension.
C
      Implicit None
      Integer iOpt
      Integer mySP(*)
      Integer N_mySP
#include "cholesky.fh"

      N_mySP = 0
      If (iOpt .eq. 1) Then
         Call Cho_P_Distrib_Vec(1,nnShl,mySP,N_mySP)
      Else
         Call Cho_P_Distrib_SP_byDim(mySP,N_mySP)
      End If

      End
      SubRoutine Cho_P_Distrib_SP_byDim(mySP,N_mySP)
C
C     Thomas Bondo Pedersen, June 2007.
C
C     Determine distribution of ShellPairs by dimension.
C
      Use Para_Info, Only: MyRank, nProcs
      use ChoArr, only: nBstSh, iSP2F
      Implicit None
      Integer mySP(*)
      Integer N_mySP
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "stdalloc.fh"

      Integer, External:: Cho_iFindSmallest

      Integer iSP, iNode, n
      Integer iAB, iA, iB

      Integer, Allocatable:: ProcDim(:)

      If (Cho_Real_Par) Then

         Call mma_allocate(ProcDim,[0,nProcs-1],Label='ProcDim')
         ProcDim(:)=0

         N_mySP = 0
         Do iSP = 1,nnShl
            iAB = iSP2F(iSP)
            Call Cho_InvPck(iAB,iA,iB,.True.)
            If (iA .eq. iB) Then
               n = nBstSh(iA)*(nBstSh(iA)+1)/2
            Else
               n = nBstSh(iA)*nBstSh(iB)
            End If
            iNode = Cho_iFindSmallest(ProcDim,SIZE(ProcDim)) - 1
            ProcDim(iNode) = ProcDim(iNode) + n
            If (iNode .eq. myRank) Then
               N_mySP = N_mySP + 1
               mySP(N_mySP) = iSP
            End If
         End Do

         Call mma_deallocate(ProcDim)

      Else

         N_mySP = nnShl
         Do iSP = 1,N_mySP
            mySP(iSP) = iSP
         End Do

      End If

      End
