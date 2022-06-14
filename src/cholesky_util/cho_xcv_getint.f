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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_XCV_GetInt(irc,ListCD,l_ListCD,ListSP,l_ListSP,
     &                          NVT,l_nVT,xInt,l_Int)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: calculate integrals (CD|J) where CD belongs to the shell
C              pairs in ListCD and J runs over all Cholesky vectors.
C
C     NOTE: this will only work if nQual and iQuAB are properly set to
C     correspond to parent diagonals. Also, mySP should be be set to a
C     trivial array [i.e. mySP(i)=i].
C
      use ChoArr, only: nDim_Batch
      Implicit None
      Integer irc
      Integer l_ListCD, l_ListSP, l_NVT, l_Int
      Integer ListCD(l_ListCD)
      Integer ListSP(l_ListSP)
      Integer NVT(l_NVT)
      Real*8  xInt(l_Int)
#include "cholesky.fh"

      Integer iSym, n
      Integer iSP, iCD

      ! Set return code
      irc=0

      ! Offsets to symmetry blocks of the integrals
      n=0
      Do iSym=1,nSym
         iOff_Col(iSym)=n
         n=n+nDim_Batch(iSym)*NVT(iSym)
      End Do

      ! Check allocation of xInt
      If (n .gt. l_Int) Then
         irc=1
         Return
      End If

      ! Calculate integrals
      Call Cho_dZero(xInt,n)
      Do iSP=1,l_ListSP
         Do iCD=1,l_ListCD
            Call Cho_MCA_CalcInt_4(xInt,n,ListCD(iCD),ListSP(iSP))
         End Do
      End Do

      End
