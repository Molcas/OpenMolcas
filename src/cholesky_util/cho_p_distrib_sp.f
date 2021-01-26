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
#include "WrkSpc.fh"

      Integer  Cho_iFindSmallest
      External Cho_iFindSmallest

      Integer iSP, iNode, n
      Integer ip_Dim, l_Dim
      Integer iAB, iA, iB

      If (Cho_Real_Par) Then

         l_Dim = nProcs
         Call GetMem('ProcDim','Allo','Inte',ip_Dim,l_Dim)
         Call iZero(iWork(ip_Dim),l_Dim)

         N_mySP = 0
         Do iSP = 1,nnShl
            iAB = iSP2F(iSP)
            Call Cho_InvPck(iAB,iA,iB,.True.)
            If (iA .eq. iB) Then
               n = nBstSh(iA)*(nBstSh(iA)+1)/2
            Else
               n = nBstSh(iA)*nBstSh(iB)
            End If
            iNode = Cho_iFindSmallest(iWork(ip_Dim),l_Dim) - 1
            iWork(ip_Dim+iNode) = iWork(ip_Dim+iNode) + n
            If (iNode .eq. myRank) Then
               N_mySP = N_mySP + 1
               mySP(N_mySP) = iSP
            End If
         End Do

         Call GetMem('ProcDim','Free','Inte',ip_Dim,l_Dim)

      Else

         N_mySP = nnShl
         Do iSP = 1,N_mySP
            mySP(iSP) = iSP
         End Do

      End If

      End
