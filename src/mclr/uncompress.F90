!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************
      SubRoutine UnCompress(ArrayIn,ArrayOut,idsym)
!
!      Uncompresses the PCG vector to a orbital rotation matrix
!
!      The redundant rotations are set to zero
!
      use Constants, only: Zero, One
      use MCLR_Data, only: nDensC, nDens2, ipMat, nB
      use input_mclr, only: nSym,TimeDep,nIsh,nOrb,nRS1,nRS2,nRS3,nBas
      Implicit None
      Real*8  ArrayIn(nDensC),ArrayOut(nDens2)
      Integer idsym

      Integer dsym
      Integer Bas(8)
      Integer IndexC, i, iSym, jSym, jT, jBas, iBas, iT, Index1, Index2
      Real*8 Fact

      indexC=0
      Fact=One
      If (idsym.lt.0) Fact=-Fact
      dsym=abs(idsym)
      ArrayOut(:)=Zero
      If (TimeDep) Then
         Do i=1,nSym
            Bas(i)=nBas(i)
         End Do
      Else
         Do i=1,nSym
            Bas(i)=nB(i)
         End Do
      End If
      Do iSym=1,nSym
       Do jSym=1,nSym
        If (iEOr(iSym-1,jSym-1)+1.eq.dSym) Then
          Do jBas=1,Bas(jSym)
          If (jBas.le.nIsh(jsym)) Then
             jT=0
          Else If (jBas.le.nIsh(jsym)+nRs1(jsym)) Then
             jT=1
          Else If (jBas.le.nIsh(jsym)+nRs1(jsym)+nRs2(jsym)) Then
             jT=2
          Else If (jBas.le.nIsh(jsym)+nRs1(jsym)+nRs2(jsym)             &
     &                               +nRs3(jsym)) Then
             jT=3
          Else
             jT=4
          End If
          Do iBas=1,nOrb(iSym)
           If (iBas.le.nIsh(isym)) Then
             iT=0
           Else If (iBas.le.nIsh(isym)+nRs1(isym)) Then
             iT=1
           Else If (iBas.le.nIsh(isym)+nRs1(isym)+nRs2(isym)) Then
             iT=2
           Else If (iBas.le.nIsh(isym)+nRs1(isym)+nRs2(isym)            &
     &                                +nRs3(isym)) Then
             iT=3
           Else
             iT=4
           End If
           If (Timedep) Then
            If (iT.ne.jT) Then
             indexC=indexc+1
             Index1=ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1
             ArrayOut(Index1)=Fact*ArrayIn(indexC)
            End If
           Else
            If (iT.gt.jT) Then
             indexC=indexc+1
             Index1=ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1
             Index2=ipMat(jSym,iSym)+(iBas-1)*nOrb(jSym)+jBas-1
             ArrayOut(Index1)=Fact*ArrayIn(indexC)
             ArrayOut(Index2)=-Fact*ArrayIn(indexC)
            End If
           End If
          End Do
         End Do
        End If
       End Do
      End Do
      End SubRoutine UnCompress
