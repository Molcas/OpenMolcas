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
* Copyright (C) 1997, Anders Bernhardsson                              *
************************************************************************
      SubRoutine UnCompress2(ArrayIn,ArrayOut,dsym)
*
*************************************************************
*
*     Change from vector to kappa   notation
*                               ip
*
*     where i is occupied and p is general
*
*     used by the preconditioner.
*
*************************************************************
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
      Integer dsym
      Real*8  ArrayIn(nDensC),ArrayOut(nDens)
      indexC=0
      Fact=1.0d0
      If (dsym.lt.0) Fact=-Fact
      dsym=abs(dsym)
      call dcopy_(nDens,0.0d0,0,ArrayOut,1)
      jT=0 ! dummy initialize
      i1=0 ! dummy initialize
      j1=0 ! dummy initialize
      Do iSym=1,nSym
       Do jSym=1,nSym
        If (iEOr(iSym-1,jSym-1)+1.eq.dSym) Then
         Do jBas=1,nB(jSym)
          If (jBas.le.nIsh(jsym)) Then
             jT=0
             i1=nIsh(isym)
             j1=nIsh(jsym)
          Else If (jBas.le.nIsh(jsym)+nRs1(jsym)) Then
             jT=1
             i1=nRs1(isym)
             j1=nRs1(jsym)
          Else If (jBas.le.nIsh(jsym)+nRs2(jsym)) Then
             jT=2
             i1=nRs2(isym)
             j1=nRs2(jsym)
          Else If (jBas.le.nIsh(jsym)+nRs3(jsym)) Then
             jT=3
             i1=nRs3(isym)
             j1=nRs3(jsym)
          End If
          Do iBas=1,nOrb(iSym)
           If (iBas.le.nIsh(isym)) Then
             iT=0
             ij=0
             ji=j1
           Else If (iBas.le.nIsh(isym)+nRs1(isym)) Then
             iT=1
             ij=0
             ji=j1
             if (it.gt.jt) then
               ij=i1
               ji=0
             end if
           Else If (iBas.le.nIsh(isym)+nRs2(isym)) Then
             iT=2
             ij=0
             ji=j1
             if (it.gt.jt) then
               ij=i1
               ji=0
             end if
           Else If (iBas.le.nIsh(isym)+nRs3(isym)) Then
             iT=3
             ij=0
             ji=j1
             if (it.gt.jt) then
               ij=i1
               ji=0
             end if
           Else
             iT=4
             ij=i1
             ji=0
           End If
           If (Timedep) Then
            If (iT.ne.jT) Then
             indexC=indexc+1
             Index1=ipMat(iSym,jSym)+(jBas-1)*norb(iSym)+iBas-1
             ArrayOut(Index1)=Fact*ArrayIn(indexC)
            End If
           Else
            If (iT.gt.jT) Then
             indexC=indexc+1
             Index1=ipMat(iSym,jSym)+
     &               (jBas-1)*norb(iSym)+iBas-ij-1
             Index2=ipMat(jSym,iSym)+
     &               (iBas-1)*norb(jSym)+jBas-ji-1
             ArrayOut(Index1)=Fact*ArrayIn(indexC)
             If (iBas.le.nB(isym))
     &        ArrayOut(Index2)=Fact*ArrayIn(indexC)
            End If
           End If
          End Do
         End Do
        End If
       End Do
      End Do
      Return
      End
