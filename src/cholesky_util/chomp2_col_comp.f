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
      SubRoutine ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Vec,nVec,Buf,lBuf,
     &                           Fac,irc)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: compute columns from a set of vectors.
C
#include "implicit.fh"
      Real*8  Col(nDim,nCol), Vec(nDim,nVec), Buf(lBuf)
      Integer iCol(nCol)

      irc = 0
      If (nDim.lt.1 .or. nCol.lt.1) Return
      If (nVec .lt. 1) Then
         If (Fac .ne. 1.0D0) Call dScal_(nDim*nCol,Fac,Col,1)
         Return
      End iF

      If (nCol.gt.1 .and. lBuf.ge.nVec) Then

         NumCol = min(lBuf/nVec,nCol)
         NumBat = (nCol - 1)/NumCol + 1

         Do iBat = 1,NumBat

            If (iBat .eq. NumBat) Then
               NumC = nCol - NumCol*(NumBat-1)
            Else
               NumC = NumCol
            End If
            iC1 = NumCol*(iBat-1) + 1

            If (lBuf .lt. NumC*nVec) Then
               irc = -1
               Return
            End If

            Call ChoMP2_Col_cp(Vec,nDim,nVec,Buf,NumC,iCol(iC1))
            Call DGEMM_('N','T',nDim,NumC,nVec,
     &                 1.0D0,Vec,nDim,Buf,NumC,
     &                 Fac,Col(1,iC1),nDim)

         End Do

      Else

         Do iC = 1,nCol
            Call dGeMV_('N',nDim,nVec,
     &                 1.0D0,Vec(1,1),nDim,Vec(iCol(iC),1),nDim,
     &                 Fac,Col(1,iC),1)
         End Do

      End If

      End
