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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine SOAdd(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,
     &                  iCmp,jCmp,iShell,jShell,
     &                  AeqB,iAO,jAO,
     &                  nComp,Label,kOper)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
      Real*8 SOInt(iBas*jBas,nSOInt), PrpInt(nPrp)
      Integer   kOper(nComp)
      Character Label*8
      Logical AeqB
*
      iRout = 130
      iPrint = nPrint(iRout)
C     Call qEnter('SOAdd')
      If (iPrint.ge.99) Then
         Call RecPrt(' In SOAdd:SOInt',' ',SOInt,iBas*jBas,nSOInt)
      End If
*
      lSO = 0
      Do 100 j1 = 0, nIrrep-1
       Do 200 i1 = 1, iCmp
        If (iAnd(IrrCmp(IndS(iShell)+i1),2**j1).eq.0) Go To 200
*
*       Scatter the SO's onto lower rectangular blocks and triangular
*       diagonal blocks.
*
        Do 300 j2 = 0, j1
         j12 = iEor(j1,j2)
         If (iAnd(lOper,2**j12).eq.0) Go To 300
         jjMx = jCmp
         If (iShell.eq.jShell .and. j1.eq.j2) jjMx = i1
         Do 400 i2 = 1, jjMx
          If (iAnd(IrrCmp(IndS(jShell)+i2),2**j2).eq.0) Go To 400
          lSO = lSO + 1
          iSO1=iAOtSO(iAO+i1,j1)
          iSO2=iAOtSO(jAO+i2,j2)
*         Write (*,*) iSO1,iAO,i1,j1,iSO2,jAO,i2,j2
*
          iPnt = iPntSO(j1,j2,lOper,nbas)
          Do 500 indAO1 = 0, iBas-1
*         Diagonal block. Store only unique elements
           jBsMax = jBas - 1
           If (j1.eq.j2 .and. iSO1.eq.iSO2) jBsMax=indAO1
           Do 600 indAO2 = 0, jBsMax
            ip = indAO2*iBas + indAO1 + 1
*
*           Move one electron integral.
*
            If (j1.eq.j2) Then
*------------Diagonal symmetry block
             If (iSO1+indAO1.ge.iSO2+indAO2) Then
              Indi = iSO1+indAO1
              Indj = iSO2+indAO2
              PrpInt(iPnt + (Indi-1)*Indi/2 + Indj) =
     &        PrpInt(iPnt + (Indi-1)*Indi/2 + Indj) + SOInt(ip,lSO)
             Else
              Indj = iSO1+indAO1
              Indi = iSO2+indAO2
              PrpInt(iPnt + (Indi-1)*Indi/2 + Indj) =
     &        PrpInt(iPnt + (Indi-1)*Indi/2 + Indj) + SOInt(ip,lSO)
             End If
            Else
*------------Off-diagonal symmetry block j1>j2
             Indi = iSO1+indAO1
             Indj = iSO2+indAO2
             nRow = nBas(j1)
*            Write (*,'(11I4,E10.5)') iSO1,iAO,i1,j1,indAO1,
*    &                                iSO2,jAO,i2,j2,indAO2,
*    &                                ip,SOInt(ip,lSO)
             PrpInt(iPnt + nRow*(Indj-1) + Indi ) =
     &       PrpInt(iPnt + nRow*(Indj-1) + Indi ) + SOInt(ip,lSO)
            End If
*
 600       Continue
 500      Continue
*
 400     Continue
 300    Continue
*
 200   Continue
 100  Continue
*
      If (iPrint.ge.99) Call GetMem(' Exit SOAdd','CHECK','ALLO',
     &                              iDum,iDum)
C     Call qExit('SOAdd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_logical(AeqB)
        Call Unused_character(Label)
        Call Unused_integer_array(kOper)
      End If
      End
