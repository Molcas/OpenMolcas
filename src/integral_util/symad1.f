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
     &                  iShll,jShll,iAO,jAO,
     &                  AOInt,iBas,jBas,nIC,iIC,
     &                  SOInt,nSOInt,nOp)
************************************************************************
*                                                                      *
* Object: to transform the one-electon matrix elements from AO basis   *
*         to SO basis.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
************************************************************************
      use Basis_Info
      use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep, Prmt
      use SOAO_Info, only: iAOtSO
      use Real_Spherical, only: iSphCr
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 AOInt(iBas*jBas,iCmp,jCmp,nIC), SOInt(iBas*jBas,nSOInt)
      Integer nOp(2)
      Integer iTwoj(0:7), jIC(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*
      iRout = 133
      iPrint = nPrint(iRout)
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
         If (iAnd(lOper,iTwoj(iIrrep)).eq.0) Cycle
         jIC(iIrrep) = iIC
         iIC = iIC + 1
 10   Continue
*
      ii = iAng*(iAng+1)*(iAng+2)/6
      jj = jAng*(jAng+1)*(jAng+2)/6
*
      lSO = 0
      Do 100 j1 = 0, nIrrep-1
         xa= DBLE(iChTbl(j1,nOp(1)))
         Do 200 i1 = 1, iCmp
            If (iAOtSO(iAO+i1,j1)<0) Cycle
            iChBs = iChBas(ii+i1)
            If (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
            pae = Prmt(iOper(nOp(1)),iChBs)
*
            Do 300 j2 = 0, nIrrep-1
               j12 = iEor(j1,j2)
*
               If (iAnd(lOper,iTwoj(j12)).eq.0) Cycle
               kIC = jIC(j12)
               xb = DBLE(iChTbl(j2,nOp(2)))
               jMx = jCmp
               If (iShell.eq.jShell .and. j1.eq.j2) jMx = i1
*
               Do 400 i2 = 1, jMx
                  If (iAOtSO(jAO+i2,j2)<0) Cycle
                  lSO = lSO + 1
                  jChBs = iChBas(jj+i2)
                  If (Shells(jShll)%Transf)
     &               jChBs = iChBas(iSphCr(jj+i2))
                  pbr = Prmt(iOper(nOp(2)),jChBs)
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
      Return
      End
