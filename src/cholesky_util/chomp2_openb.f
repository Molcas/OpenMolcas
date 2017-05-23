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
      SubRoutine ChoMP2_OpenB(iOpt,iSym,iBatch)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: open (iOpt=1), close and keep (iOpt=2), or close and
C              delete (iOpt=3) Cholesky vector files for MP2 program
C              (batch vectors).
C              For iOpt=0, the units are initialized (to -1).
C
#include "implicit.fh"
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'ChoMP2_OpenB')

      Character*2 BaseNm
      Parameter (BaseNm = '_I')

      Character*6 BtchNm

      LnT1am(i,j)=iWork(ip_LnT1am-1+nSym*(j-1)+i)
      lUnit(i,j)=iWork(ip_lUnit-1+nSym*(j-1)+i)

C     Initialize units and return for iOpt=0.
C     ---------------------------------------

      If (iOpt .eq. 0) Then
         iWork(ip_lUnit-1+nSym*(iBatch-1)+iSym) = -1
         Return
      End If

C     Open or close files.
C     --------------------

      If (iOpt .eq. 1) Then
         If (LnT1am(iSym,iBatch) .gt. 0) Then
            If (iBatch .lt. 10) Then
               Write(BtchNm,'(A2,I1,A2,I1)') BaseNm,iSym,'__',iBatch
            Else If (iBatch .lt. 100) Then
               Write(BtchNm,'(A2,I1,A1,I2)') BaseNm,iSym,'_',iBatch
            Else If (iBatch .lt. 1000) Then
               Write(BtchNm,'(A2,I1,I3)')    BaseNm,iSym,iBatch
            Else ! note: due to restriction in filename length...
               Call ChoMP2_Quit(SecNam,'Too many batches',
     &                          '(Current max. is 999)')
               BtchNm = '?!?!?!' ! too avoid compiler warnings...
            End If
            lU = 7
            Call daName_MF_WA(lU,BtchNm)
         Else
            lU = -1
         End If
         iWork(ip_lUnit-1+nSym*(iBatch-1)+iSym) = lU
      Else If (iOpt .eq. 2) Then
         lU = lUnit(iSym,iBatch)
         If (lU .gt. 0) Then
            Call daClos(lU)
            iWork(ip_lUnit-1+nSym*(iBatch-1)+iSym) = -1
         End If
      Else If (iOpt .eq. 3) Then
         lU = lUnit(iSym,iBatch)
         If (lU .gt. 0) Then
            Call daEras(lU)
            iWork(ip_lUnit-1+nSym*(iBatch-1)+iSym) = -1
         End If
      Else
         Call ChoMP2_Quit(SecNam,'iOpt out of bounds',' ')
      End If

      End
