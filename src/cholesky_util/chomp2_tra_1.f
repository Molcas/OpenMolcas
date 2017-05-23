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
      SubRoutine ChoMP2_Tra_1(COcc,CVir,Diag,DoDiag,Wrk,lWrk,iSym)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: transform Cholesky vectors to (ai) MO basis for symmetry
C              block iSym. Files are assumed open.
C              If requested (DoDiag=.true.), compute (ai|ai) integral
C              diagonal.
C
#include "implicit.fh"
      Real*8  COcc(*), CVir(*), Diag(*), Wrk(lWrk)
      Logical DoDiag
#include "cholesky.fh"
#include "choptr.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'ChoMP2_Tra_1')

      Integer  Cho_lRead
      External Cho_lRead

      Integer ai

      Parameter (N2 = InfVec_N2)

      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)

C     Check if anything to do.
C     ------------------------

      If (NumCho(iSym).lt.1 .or. nT1am(iSym).lt.1) Return

C     Initialize Diag (if needed).
C     ----------------------------

      If (DoDiag) Call Cho_dZero(Diag,nT1am(iSym))

C     Allocate memory for half-transformed vector.
C     --------------------------------------------

      lHlfTr = nT1AOT(iSym)

      kHlfTr = 1
      kEnd0  = kHlfTr + lHlfTr
      lWrk0  = lWrk   - kEnd0 + 1
      If (lWrk0 .lt. (nT1am(iSym)+nnBstR(iSym,1))) Then
         Call ChoMP2_Quit(SecNam,'insufficient memory','[0]')
      End If

C     Reserve memory for reading AO vectors.
C     --------------------------------------

      lRead = Cho_lRead(iSym,lWrk0)
      If (lRead .lt. 1) Then
         Write(6,*) SecNam,': memory error: lRead = ',lRead
         Call ChoMP2_Quit(SecNam,'memory error',' ')
         lWrk1 = 0 ! to avoid compiler warnings...
      Else
         lWrk1 = lWrk0 - lRead
         If (lWrk1 .lt. nT1am(iSym)) Then
            lWrk1 = nT1am(iSym)
            lRead = lWrk - nT1am(iSym)
         End If
      End If

C     Set up batch.
C     -------------

      nMOVec = min(lWrk1/nT1am(iSym),NumCho(iSym))
      If (nMOVec .lt. 1) Then
         Call ChoMP2_Quit(SecNam,'insufficient memory','[1]')
      End If
      NumBat = (NumCho(iSym) - 1)/nMOVec + 1

C     Set reduced set handles.
C     ------------------------

      iRedC = -1
      iLoc  = 3

C     Transform each batch of vectors and compute diagonal contributions
C     (if requested).
C     ------------------------------------------------------------------

      Do iBat = 1,NumBat

         If (iBat .eq. NumBat) Then
            NumV = NumCho(iSym) - nMOVec*(NumBat - 1)
         Else
            NumV = nMOVec
         End If
         iVec1 = nMOVec*(iBat - 1) + 1
         iVec2 = iVec1 + NumV - 1

         lChoMO = nT1am(iSym)*NumV

         kChoMO = kEnd0
         kChoAO = kChoMO + lChoMO
         lChoAO = lWrk0  - kChoAO + 1

         kOffMO = kChoMO
         jVec1  = iVec1
         Do While (jVec1 .le. iVec2)

            jNum = 0
            Call Cho_VecRd(Wrk(kChoAO),lChoAO,jVec1,iVec2,iSym,
     &                     jNum,iRedC,mUsed)
            If (jNum .lt. 1) Then
               Call ChoMP2_Quit(SecNam,
     &                          'insufficient memory','[2]')
            End If

            kOff = kChoAO
            Do jVec = 1,jNum
               iVec = jVec1 + jVec - 1
               iRed = InfVec(iVec,2,iSym)
               If (iRedC .ne. iRed) Then
                  irc = 0
                  Call Cho_X_SetRed(irc,iLoc,iRed)
                  If (irc .ne. 0) Then
                     Call ChoMP2_Quit(SecNam,'error in Cho_X_SetRed',
     &                                ' ')
                  End If
                  iRedC = iRed
               End If
               Call ChoMP2_TraVec(Wrk(kOff),Wrk(kOffMO),COcc,CVir,
     &                            Wrk(kHlfTr),lHlfTr,iSym,1,1,iLoc)
               kOff   = kOff   + nnBstR(iSym,iLoc)
               kOffMO = kOffMO + nT1am(iSym)
            End Do

            jVec1 = jVec1 + jNum

         End Do

         iOpt = 1
         iAdr = nT1am(iSym)*(iVec1 - 1) + 1
         Call ddaFile(lUnit_F(iSym,1),iOpt,Wrk(kChoMO),lChoMO,iAdr)

         If (DoDiag) Then
            Do iVec = 1,NumV
               kOff = kChoMO + nT1am(iSym)*(iVec-1) - 1
               Do ai = 1,nT1am(iSym)
                  Diag(ai) = Diag(ai) + Wrk(kOff+ai)*Wrk(kOff+ai)
               End Do
            End Do
         End If

      End Do

      End
