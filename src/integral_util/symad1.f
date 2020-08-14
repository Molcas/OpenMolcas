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
      Subroutine SymAd1(lOper,iAng,jAng,iCmp,jCmp,iShell,jShell,
     &                  iShll,jShll,AOInt,iBas,jBas,nIC,iIC,
     &                  SOInt,nSOInt,nOp)
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
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
      Real*8 AOInt(iBas*jBas,iCmp,jCmp,nIC), SOInt(iBas*jBas,nSOInt)
      Integer nOp(2)
      Real*8 Prmt(0:7)
      Integer iTwoj(0:7), jIC(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement functions
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iRout = 133
      iPrint = nPrint(iRout)
*     Call qEnter('SymAd1')
      If (iPrint.ge.99) Then
         Write (6,*) ' lOper=',lOper
         Write (6,*) ' nSOInt=',nSOInt
         Call RecPrt(' In SymAd1: AOInt',' ',AOInt,iBas*jBas,
     &                iCmp*jCmp*nIC)
         Call RecPrt(' In SymAd1: SOInt',' ',SOInt,iBas*jBas,nSOInt)
         Write (6,*) ' iIC=',iIC
      End If
      Do 10 iIrrep = 0, nIrrep-1
         jIC(iIrrep) = -999999999
         If (iAnd(lOper,iTwoj(iIrrep)).eq.0) Go To 10
         jIC(iIrrep) = iIC
         iIC = iIC + 1
 10   Continue
*
      ii = iAng*(iAng+1)*(iAng+2)/6
      jj = jAng*(jAng+1)*(jAng+2)/6
*
      lSO = 0
      Do 100 j1 = 0, nIrrep-1
         xa= rChTbl(j1,nOp(1))
         Do 200 i1 = 1, iCmp
            If (iAnd(IrrCmp(IndS(iShell)+i1),iTwoj(j1)).eq.0) Go To 200
            iChBs = iChBas(ii+i1)
            If (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
            pae = xPrmt(iOper(nOp(1)),iChBs)
*
*old        Do 300 j2 = 0, j1
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
                  jChBs = iChBas(jj+i2)
                  If (Shells(jShll)%Transf)
     &               jChBs = iChBas(iSphCr(jj+i2))
                  pbr = xPrmt(iOper(nOp(2)),jChBs)
                  Call DaXpY_(iBas*jBas,xa*pae*xb*pbr,
     &                       AOInt(1,i1,i2,kIC),1,
     &                       SOInt(1,lSO),1)
 400           Continue
*
 300        Continue
*
 200     Continue
 100  Continue
      If (lSO.ne.nSOInt) Then
         Call WarningMessage(2,'Error in SymAd1, lSO.ne.nSOInt')
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
      End
