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
* Copyright (C) 2000, Roland Lindh                                     *
*               Ajitha Devarajan                                       *
************************************************************************
      Subroutine DFT_IntX(Do_NInt_d,Do_NInt,
     &                    Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                    FckInt,nFckInt,SOTemp,nSOTemp,
     &                    TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                    nSym,iSpin,Flop,Rho,nRho,Scr,nScr,
     &                    Fact,ndc,mAO,TabAOMax,T_X,list_bas,nFn)
************************************************************************
*                                                                      *
* Object: to compute contributions to                                  *
*                                                                      *
*         <m|dF/drho|n> ; integrals over the potential                 *
*                                                                      *
*         where                                                        *
*                                                                      *
*         F(r)=rho(r)*e(rho(r),grad[rho(r)])                           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
*             D.Ajitha:Modifying for the new Kernel outputs            *
************************************************************************
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
      External Do_NInt_d, Do_NInt
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 Weights(mGrid), SOTemp(nSOTemp,iSpin), Fact(ndc**2),
     &       TabAO(nTabAO), Scr(nScr*mGrid), Rho(nRho,mGrid),
     &       AOInt(nAOInt*nAOInt,iSpin), FckInt(nFckInt,iSpin),
     &       dF_dRho(ndF_dRho,mGrid), TabAOMax(nlist_s)
      Integer nOp(2), list_s(2,nlist_s), ipTabAO(nlist_s),
     &        list_bas(2,nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
C     Call QEnter('DFT_IntX')
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the desired AO integrand here from the AOs, accumulate
*     contributions to the SO integrals on the fly.
*
*define _DEBUG_
#ifdef _DEBUG_
      Debug=.True.
#endif
      VMax=0.0D0
      iSmLbl=1
*
      nGrid_Tot=0
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         TMax_i=TabAOMax(ilist_s)
         If (TMax_i.le.T_X) Go To 999
         kDCRE = list_s(2,ilist_s)
         iShll = iSD( 0,iSkal)
         iAng  = iSD( 1,iSkal)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         iAO   = iSD( 7,iSkal)
         IndShl= iSD( 8,iSkal)
         mdci  = iSD(10,iSkal)
         iShell= iSD(11,iSkal)
*
         nOp(1) = NrOpr(kDCRE,iOper,nIrrep)
*
         Do jlist_s=ilist_s,nlist_s
            jSkal = list_s(1,jlist_s)
            TMax_j=TabAOMax(jlist_s)
            If (TMax_j.le.T_X) Go To 998
            kDCRR = list_s(2,jlist_s)
            jShll = iSD( 0,jSkal)
            jAng  = iSD( 1,jSkal)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            jAO   = iSD( 7,jSkal)
            JndShl= iSD( 8,jSkal)
            mdcj  = iSD(10,jSkal)
            jShell= iSD(11,jSkal)
*
            nOp(2) = NrOpr(kDCRR,iOper,nIrrep)
*
            nSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,IndShl,JndShl)
*
            ij = (mdcj-1)*ndc + mdci
*                                                                      *
************************************************************************
*                                                                      *
            If (ilist_s.eq.jlist_s) Then
               Call Do_NInt_d(AOInt,nAOInt,ndF_dRho, dF_dRho,
     &                        Weights,mGrid,Rho,nRho,
     &                        Scr, TabAO(ipTabAO(iList_s)),iCmp,
     &                        iBas_Eff,nGrid_Tot,iSpin,mAO,nFn)
               jx = iDAMax_(iSpin*mAO*mGrid*iBas_Eff*iCmp,Scr,1)
               VMax=Abs(Scr(jx))
            Else
               If (VMax*TMax_j.le.T_X) Go To 998
               Call Do_NInt(AOInt,nAOInt,mGrid,
     &                      Scr,iCmp,iBas_Eff,
     &                      TabAO(ipTabAO(jList_s)),jCmp,jBas_Eff,
     &                      nGrid_Tot,iSpin,mAO,nFn)
            End If
*
            Do iD = 1, iSpin
*
#ifdef _DEBUG_
               If (Debug) Then
                  nAOInt_j=jBas_Eff*jCmp
                  nAOInt_i=iBas_Eff*iCmp
                  mAOInt=nAOInt_i*nAOInt_j
                  Call RecPrt('Kernel: AOInt',' ',AOInt(1,iD),1,mAOInt)
               End If
#endif
               If (nIrrep.ne.1) Then
*
*---------------- Distribute contributions to the SO integrals
*
                  nIC = 1
                  iIC = 1
                  Call FZero(SOTemp(1,iD),iBas*jBas*nSO)
                  Call SymAdd2(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                         iShell,jShell,iShll,jShll,
     &                         IndShl,JndShl,AOInt(1,iD),
     &                         iBas,iBas_Eff,jBas,jBas_Eff,
     &                         nIC,iIC,SOTemp(1,iD),nSO,nOp,
     &                         iSkal,jSkal)
*
*---------------- Here scatter the result
*
                  If (Fact(ij).ne.One)
     &               Call DScal_(nSO*iBas*jBas,Fact(ij),SOTemp(1,iD),1)
*
                  Call SOAdd(SOTemp(1,iD),iBas,jBas,nSO,
     &                       FckInt(1,iD),nFckInt,iSmLbl,
     &                       iCmp,jCmp,iShell,jShell,IndShl,JndShl,
     &                       iSkal.eq.jSkal,iAO,jAO)
*
               Else
*
                  Call AOAdd_nq(AOInt(1,iD),iBas,iBas_Eff,jBas,jBas_Eff,
     &                          FckInt(1,iD),nFckInt,
     &                          iCmp,jCmp,iShell,jShell,iAO,jAO)

*
               End If
*
            End Do                  ! iD
*                                                                      *
************************************************************************
*                                                                      *
*
 998        Continue
         End Do                     ! jlist_s
 999     Continue
      End Do                        ! ilist_s
      Flop=Flop+DBLE(nGrid_Tot)
*
C     Call QExit('DFT_Int1')
#ifdef _DEBUG_
      Debug=.False.
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nSym)
      End
