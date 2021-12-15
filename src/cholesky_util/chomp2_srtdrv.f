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
      SubRoutine ChoMP2_SrtDrv(irc,DelOrig)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: presort Cholesky vectors according to batch structure in
C              MP2 program.
C
C     DelOrig: input : flag for deleting original vector files.
C              output: flag to tell that at least 1 symmetry block has
C                      in fact been deleted.
C
#include "implicit.fh"
      Logical DelOrig
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character*6  ThisNm
      Character*13 SecNam
      Parameter (SecNam = 'ChoMP2_SrtDrv', ThisNm = 'SrtDrv')

      LnT1am(i,j)=iWork(ip_LnT1am-1+nSym*(j-1)+i)
      lUnit(i,j)=iWork(ip_lUnit-1+nSym*(j-1)+i)

      irc = 0
      If (nBatch .lt. 1) Return


C     Allocate available memory.
C     --------------------------

      Call GetMem('GetMax','Max ','Real',ipWrk,lWrk)
      Call GetMem('SrtMax','Allo','Real',ipWrk,lWrk)

C     Set vector type (i.e., transformed vectors or vectors from (ai|bj)
C     decomposition. Decide whether original files should be deleted.
C     ------------------------------------------------------------------

      If (DecoMP2) Then
         iTyp = 2
      Else
         iTyp = 1
      End If

      If (DelOrig) Then
         iClos = 3
      Else
         iClos = 2
      End If
      DelOrig = .false.

C     Start symmetry loop.
C     --------------------

      Do iSym = 1,nSym

C        Set number of vectors.
C        ----------------------

         If (iTyp .eq. 1) Then
            nSrtVec = NumCho(iSym)
         Else If (iTyp .eq. 2) Then
            nSrtVec = nMP2Vec(iSym)
         Else
            irc = -1
            Go To 1 ! exit
         End If

         If (nT1am(iSym).gt.0 .and. nSrtVec.gt.0) Then

C           Set up vector batch.
C           --------------------

            LnT1amx = 0
            Do iBatch = 1,nBatch
               LnT1amx = max(LnT1amx,LnT1am(iSym,iBatch))
            End Do

            MinMem = nT1am(iSym) + LnT1amx
            NumVec = min(lWrk/MinMem,nSrtVec)
            If (NumVec .lt. 1) Then
               irc = 1
               Go To 1 ! exit
            Else
               nBat = (nSrtVec - 1)/NumVec + 1
            End If

C           Open full vector file.
C           ----------------------

            Call ChoMP2_OpenF(1,iTyp,iSym)

C           Start batch loop.
C           -----------------

            Do iBat = 1,nBat

               If (iBat .eq. nBat) Then
                  NumV = nSrtVec - NumVec*(nBat-1)
               Else
                  NumV = NumVec
               End If
               iVec1 = NumVec*(iBat-1) + 1

C              Read batch of vectors.
C              ----------------------

               lChoMO = nT1am(iSym)*NumV

               kChoMO = ipWrk

               iOpt = 2
               lTot = lChoMO
               iAdr = nT1am(iSym)*(iVec1-1) + 1
               Call ddaFile(lUnit_F(iSym,iTyp),iOpt,Work(kChoMO),lChoMO,
     &                      iAdr)

C              Sort and write to disk.
C              -----------------------

               kSort  = ipWrk + lChoMO
               lSort  = lWrk  - lChoMO
               Do iBatch = 1,nBatch
                  lTot = LnT1am(iSym,iBatch)*NumV
                  If (lSort .lt. lTot) Then
                     Call ChoMP2_Quit(SecNam,'sort batch error','[0]')
                  End If
                  Call ChoMP2_Srt(Work(kChoMO),Work(kSort),NumV,
     &                            iSym,iBatch)
                  Call ChoMP2_OpenB(1,iSym,iBatch)
                  iOpt = 1
                  iAdr = LnT1am(iSym,iBatch)*(iVec1-1) + 1
                  Call ddaFile(lUnit(iSym,iBatch),iOpt,Work(kSort),lTot,
     &                         iAdr)
                  Call ChoMP2_OpenB(2,iSym,iBatch)
               End Do

            End Do

C           Close (and possibly delete) full vector file.
C           ---------------------------------------------

            Call ChoMP2_OpenF(iClos,iTyp,iSym)
            DelOrig = iClos .eq. 3

         End If

      End Do

    1 Call GetMem('SrtMax','Free','Real',ipWrk,lWrk)
      End
