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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine Compress(ArrayIn,ArrayOut,dsym)
*
*      Compresses the orbital rotation matrix to
*      the vector that is used in the PCG routines
*      the indexes are ordered to fit the preconditioner
*
*      The redundant rotations are set to zero
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
      Integer dsym
      Real*8  ArrayIn(nDens),ArrayOut(nDensC)
      indexC=0
      call dcopy_(nDensC,[0.0d0],0,ArrayOut,1)
      Do iSym=1,nSym
       Do jSym=1,nSym
        If (iEOr(iSym-1,jSym-1)+1.eq.abs(dSym)) Then
         Do jBas=1,nOrb(jSym)
          If (jBas.le.nIsh(jsym)) Then
             jT=0
          Else If (jBas.le.nIsh(jsym)+nRs1(jsym)) Then
             jT=1
          Else If (jBas.le.nIsh(jsym)+nRs2(jsym)) Then
             jT=2
          Else If (jBas.le.nIsh(jsym)+nRs3(jsym)) Then
             jT=3
          Else
             jT=4
          End If
          Do iBas=1,nOrb(iSym)
           If (iBas.le.nIsh(isym)) Then
             iT=0
           Else If (iBas.le.nIsh(isym)+nRs1(isym)) Then
             iT=1
           Else If (iBas.le.nIsh(isym)+nRs2(isym)) Then
             iT=2
           Else If (iBas.le.nIsh(isym)+nRs3(isym)) Then
             iT=3
           Else
             iT=4
           End If
           If (TimeDep) Then
            If (iT.ne.jT) Then !
             indexC=indexc+1 !
             Index1=ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1 !
             Index2=ipMat(jSym,iSym)+(iBas-1)*nOrb(jSym)+jBas-1 !
             ArrayOut(IndexC)=ArrayIn(index1) !
            End If !
           Else
            If (iT.gt.jT) Then
             indexC=indexc+1
             Index1=ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1
             ArrayOut(IndexC)=ArrayIn(index1)
            End If
           End If
          End Do
         End Do
        End If
       End Do
      End Do
      If (indexc.ne.ndensc) Call SysAbendMsg('compress',
     & 'indexc.ne.ndensc',' ')
      Return
      End
