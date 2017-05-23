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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_DecChk_Int(irc,lUnit,Col,Nai,Nbj,ibj1,NumVec,
     &                             Work,lWork,Fac)
C
C     Thomas Bondo Pedersen, Jan. 2005.
C
C     Purpose: compute consecutive columns of (ai|bj) matrix from
C              vectors on file (unit: lUnit)
C              vectors from MP2 decomposition).
C
#include "implicit.fh"
      Real*8  Col(Nai,Nbj),Work(lWork)

      irc = 0

C     Check dimensions.
C     -----------------

      If (Nai.lt.1 .or. Nbj.lt.1 .or. Nbj.gt.Nai) Then
         irc = -1
         return
      End If
      ibj2 = ibj1 + Nbj - 1
      If (ibj1.lt.1 .or. ibj2.gt.Nai) Then
         irc = -2
         Return
      End If

C     Scale result array.
C     -------------------

      Call dScal_(Nai*Nbj,Fac,Col,1)
      If (NumVec .lt. 1) Return

C     Set up batch.
C     -------------
      nVec = min(lWork/Nai,NumVec)
      If (nVec .lt. 1) Then
         irc = 1
         Return
      End If
      nBat = (NumVec - 1)/nVec + 1

C     Start batch loop.
C     -----------------

      Do iBat = 1,nBat

C        Set batch info.
C        ---------------

         If (iBat .eq. nBat) Then
            NumV = NumVec - nVec*(nBat - 1)
         Else
            NumV = nVec
         End If
         iVec1 = nVec*(iBat - 1) + 1

C        Read vectors.
C        -------------

         iOpt = 2
         lTot = Nai*NumV
         iAdr = Nai*(iVec1 - 1) + 1
         Call ddaFile(lUnit,iOpt,Work,lTot,iAdr)

C        Compute integrals.
C        ------------------

         Call DGEMM_('N','T',Nai,Nbj,NumV,
     &              1.0d0,Work,Nai,Work(ibj1),Nai,
     &              1.0d0,Col,Nai)

      End Do

      End
