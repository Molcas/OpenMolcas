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
      Subroutine SymAdd2(lOper,iAng,jAng,iCmp,jCmp,iShell,jShell,
     &                   iShll,jShll,AOInt,iBas,iBas_Eff,
     &                                     jBas,jBas_Eff,nIC,iIC,
     &                   SOInt,nSOInt,nOp,iSkal,jSkal)
************************************************************************
*                                                                      *
* Object: to transform the one-electon matrix elements from AO basis   *
*         to SO basis.                                                 *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DaXpY   (ESSL)                                          *
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
      Real*8 AOInt(iBas_Eff*jBas_Eff,iCmp,jCmp,nIC),
     &       SOInt(iBas*jBas,nSOInt)
      Integer nOp(2)
      Real*8 Prmt(0:7)
      Integer iTwoj(0:7), jIC(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
      iRout = 133
      iPrint = nPrint(iRout)
*     Call qEnter('SymAd1')
      If (iPrint.ge.99) Then
         Write (6,*) ' lOper=',lOper
         Write (6,*) ' nSOInt=',nSOInt
         Call RecPrt(' In SymAdd: AOInt',' ',AOInt,iBas*jBas,
     &                iCmp*jCmp*nIC)
         Call RecPrt(' In SymAdd: SOInt',' ',SOInt,iBas*jBas,nSOInt)
         Write (6,*) ' iIC=',iIC
      End If
      Do 10 iIrrep = 0, nIrrep-1
         jIC(iIrrep) = -999999999
         If (iAnd(lOper,iTwoj(iIrrep)).eq.0) Go To 10
         jIC(iIrrep) = iIC
         iIC = iIC + 1
 10   Continue
*
      lSO = 0
      iAdd = iBas-iBas_Eff
      jAdd = jBas-jBas_Eff
      Do 100 j1 = 0, nIrrep-1
         xa = rChTbl(j1,nOp(1))
         Do 200 i1 = 1, iCmp
            If (iAnd(IrrCmp(IndS(iShell)+i1),iTwoj(j1)).eq.0) Go To 200
*
            Do 300 j2 = 0, nIrrep-1
               j12 = iEor(j1,j2)
*
               If (iAnd(lOper,iTwoj(j12)).eq.0) Go To 300
               kIC = jIC(j12)
               xb = rChTbl(j2,nOp(2))
               jMx = jCmp
               If (iShell.eq.jShell .and. j1.eq.j2) jMx = i1
*
               Do 400 i2 = 1, jMx
                  If (iAnd(IrrCmp(IndS(jShell)+i2),iTwoj(j2)).eq.0)
     &               Go To 400
                  lSO = lSO + 1
*
                  Do iB_Eff = 1, iBas_Eff
                     iB = iB_Eff + iAdd
                     Do jB_Eff = 1, jBas_Eff
                        jB = jB_Eff + jAdd
*
                        iFrom=(jB_Eff-1)*iBas_Eff+iB_Eff
                        iTo  =(jB    -1)*iBas    +iB
                        SOInt(iTo,lSO)=SOInt(iTo,lSO)
     &                                +xa*xb*AOInt(iFrom,i1,i2,kIC)
                        If (iSkal.eq.jSkal.and.nOp(1).ne.nOp(2)) Then
                           iTo  =(iB    -1)*jBas    +jB
                           SOInt(iTo,lSO)=SOInt(iTo,lSO)
     &                                   +xa*xb*AOInt(iFrom,i2,i1,kIC)
                        End If
*
                     End Do
                  End Do
*
 400           Continue
*
 300        Continue
*
 200     Continue
 100  Continue
      If (lSO.ne.nSOInt) Then
         Call WarningMessage(2,'Error in SymAdd, lSO.ne.nSOInt')
         Call Abend()
      End If
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In SymAd1: SOInt',' ',SOInt,iBas*jBas,nSOInt)
      End If
      If (iPrint.ge.59) Call GetMem(' Exit SymAd1','CHECK','REAL',
     &                              iDum,iDum)
*     Call qExit('SymAd1')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iAng)
         Call Unused_integer(jAng)
         Call Unused_integer(iShll)
         Call Unused_integer(jShll)
      End If
      End
