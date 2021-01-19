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
      SubRoutine ChoMP2g_Tra_1(COrb1,COrb2,Diag,DoDiag,Wrk,lWrk,iSym,
     &                         iMoType1, iMoType2)
C
C     Thomas Bondo Pedersen, Dec. 2010.
C
C     Purpose: transform Cholesky vectors to (pq) MO basis for symmetry
C              block iSym. Files are assumed open.
C              If requested (DoDiag=.true.), compute (pq|pq) integral
C              diagonal.
C
      use ChoSwp, only: InfVec
#include "implicit.fh"
      Real*8  COrb1(*), COrb2(*), Diag(*), Wrk(lWrk)
      Logical DoDiag
#include "cholesky.fh"
#include "choptr.fh"
#include "chomp2.fh"
#include "chomp2g.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'ChoMP2_Tra_1')

      Integer  Cho_lRead
      External Cho_lRead

      Integer pq

*     Check what type of Cholesky vector to make (fro-occ, occ-occ.....)
      iVecType = iMoType2 + (iMoType1-1)*nMoType

C     Check if anything to do.
C     ------------------------

      If (NumCho(iSym).lt.1 .or. nMoMo(iSym,iVecType).lt.1) Return

C     Initialize Diag (if needed).
C     ----------------------------

      If (DoDiag) Call Cho_dZero(Diag,nMoMo(iSym,iVecType))

C     Allocate memory for half-transformed vector.
C     --------------------------------------------

      lHlfTr = nMoAo(iSym,iMoType1)

      kHlfTr = 1
      kEnd0  = kHlfTr + lHlfTr
      lWrk0  = lWrk   - kEnd0 + 1
      If (lWrk0 .lt. (nMoMo(iSym,iVecType)+nnBstR(iSym,1))) Then
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
         If (lWrk1 .lt. nMoMo(iSym,iVecType)) Then
            lWrk1 = nMoMo(iSym,iVecType)
            lRead = lWrk - nMoMo(iSym,iVecType)
         End If
      End If

C     Set up batch.
C     -------------

      nMOVec = min(lWrk1/nMoMo(iSym,iVecType),NumCho(iSym))
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

         lChoMO = nMoMo(iSym,iVecType)*NumV

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
*
               Call ChoMP2g_TraVec(Wrk(kOff),Wrk(kOffMO),COrb1,COrb2,
     &                            Wrk(kHlfTr),lHlfTr,iSym,1,1,iLoc,
     &                            iMoType1,iMoType2)
               kOff   = kOff   + nnBstR(iSym,iLoc)
               kOffMO = kOffMO + nMoMo(iSym,iVecType)
            End Do

            jVec1 = jVec1 + jNum

         End Do


         iOpt = 1
         iAdr = nAdrOff(iSym) +
     &          nMoMo(iSym,iVecType)*(iVec1 - 1) + 1
         iAdrOff(iSym,iVecType) = nAdrOff(iSym)
         Call ddaFile(lUnit_F(iSym,1),iOpt,Wrk(kChoMO),lChoMO,iAdr)

         If (DoDiag) Then
            Do iVec = 1,NumV
               kOff = kChoMO + nMoMo(iSym,iVecType)*(iVec-1) - 1
               Do pq = 1,nMoMo(iSym,iVecType)
                  Diag(pq) = Diag(pq) + Wrk(kOff+pq)*Wrk(kOff+pq)
               End Do
            End Do
         End If

      End Do
*     When we reach this point we have written all vectors of the present
*     type for this symmetry and need to remember were we should continue
*     to write the next type.
      If(iVecType.ne.9) Then
         nAdrOff(iSym) = iAdr-1
      End If
      End
