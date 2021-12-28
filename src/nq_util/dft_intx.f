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
*#define _NEWCODE_
#ifdef _NEWCODE_
      Subroutine DFT_IntX(Do_NInt_d,Do_NInt
     &                    Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                    FckInt,nFckInt,
     &                    ipTabAO,dF_dRho,ndF_dRho,
     &                    nSym,iSpin,Flop,Scr,nScr,
     &                    Fact,ndc,mAO,list_bas,nFn)
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
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
*             D.Ajitha:Modifying for the new Kernel outputs            *
************************************************************************
      use iSD_data
      use Symmetry_Info, only: nIrrep
      use nq_Grid, only: TabAO_Pack, TabAO, Grid_AO, iBfn_Index,
     &                   AOIntegrals => Dens_AO
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (A-H,O-Z)
      External Do_NInt_d, Do_NIntX
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 Weights(mGrid), Fact(ndc**2), Scr(nScr*mGrid),
     &       AOInt(nAOInt*nAOInt,iSpin), FckInt(nFckInt,iSpin),
     &       dF_dRho(ndF_dRho,mGrid)
      Integer nOp(2), list_s(2,nlist_s), ipTabAO(nlist_s),
     &        list_bas(2,nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the desired AO integrand here from the AOs, accumulate
*     contributions to the SO integrals on the fly.
*
      nGrid_Tot=0
*
      nBfn = Size(Grid_AO,3)
      Call Do_NInt_d(ndF_dRho, dF_dRho,Weights,mGrid,
     &               Grid_AO, TabAO,nBfn,nGrid_Tot,iSpin,mAO,nFn)
      Call Do_NIntX(AOIntegrals,mGrid,Grid_AO,TabAO,nBfn,
     &                   nGrid_Tot,iSpin,mAO,nFn)
*                                                                      *
************************************************************************
*                                                                      *
*     Set up an indexation translation between the running index of
*     the AOIntegrals and the actuall basis function index
*
      iBfn = 0
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         kDCRE = list_s(2,ilist_s)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         iAO   = iSD( 7,iSkal)
         mdci  = iSD(10,iSkal)

         iAdd = iBas-iBas_Eff
         Do i1 = 1, iCmp
            iSO1 = iAOtSO(iAO+i1,0)
            Do i2 = 0, iBas_Eff-1
               IndAO1 = i2 + iAdd
               Indi = iSO1 + indAO1

               iBfn = iBfn + 1
               iBfn_Index(1,iBfn) = Indi
               iBfn_Index(2,iBfn) = ilist_s
               iBfn_Index(3,iBfn) = i1
               iBfn_Index(4,iBfn) = i2+1
               iBfn_Index(5,iBfn) = mdci
               iBfn_Index(6,iBfn) = IndAO1
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
     Do iD = 1, iSpin
     If (nIrrep.eq.1) Then
        Call AOAdd_Full(AOIntegrals(:,:,iD),nBfn,PrpInt(:,iD),nPrp)
     Else
     End If
     End Do
*                                                                      *
************************************************************************
*                                                                      *
      Flop=Flop+DBLE(nGrid_Tot)
*
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(nSym)
#endif
      End
#else
      Subroutine DFT_IntX(Do_NInt_d,Do_NInt,
     &                    Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                    FckInt,nFckInt,
     &                    ipTabAO,dF_dRho,ndF_dRho,
     &                    nSym,iSpin,Flop,Scr,nScr,
     &                    Fact,ndc,mAO,list_bas,nFn)
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
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
*             D.Ajitha:Modifying for the new Kernel outputs            *
************************************************************************
      use iSD_data
      use Symmetry_Info, only: nIrrep
      use nq_Grid, only: TabAO_Pack
      Implicit Real*8 (A-H,O-Z)
      External Do_NInt_d, Do_NInt
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 Weights(mGrid), Fact(ndc**2), Scr(nScr*mGrid),
     &       AOInt(nAOInt*nAOInt,iSpin), FckInt(nFckInt,iSpin),
     &       dF_dRho(ndF_dRho,mGrid)
      Integer nOp(2), list_s(2,nlist_s), ipTabAO(nlist_s),
     &        list_bas(2,nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the desired AO integrand here from the AOs, accumulate
*     contributions to the SO integrals on the fly.
*
      nGrid_Tot=0
*
      iSmLbl=1
*                                                                      *
************************************************************************
*                                                                      *
*     The purpose of this loop is to symmetry adapt the contributions.
*                                                                      *
************************************************************************
*                                                                      *
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         kDCRE = list_s(2,ilist_s)
         iShll = iSD( 0,iSkal)
         iAng  = iSD( 1,iSkal)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         iAO   = iSD( 7,iSkal)
         mdci  = iSD(10,iSkal)
         iShell= iSD(11,iSkal)
*
         nOp(1) = NrOpr(kDCRE)
         Call Do_NInt_d(ndF_dRho, dF_dRho,
     &                  Weights,mGrid,
     &                  Scr, TabAO_Pack(ipTabAO(iList_s):),
     &                  iCmp*iBas_Eff,nGrid_Tot,iSpin,mAO,nFn)
*
         Do jlist_s=ilist_s,nlist_s
            jSkal = list_s(1,jlist_s)
            kDCRR = list_s(2,jlist_s)
            jShll = iSD( 0,jSkal)
            jAng  = iSD( 1,jSkal)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            jAO   = iSD( 7,jSkal)
            mdcj  = iSD(10,jSkal)
            jShell= iSD(11,jSkal)
*
            nOp(2) = NrOpr(kDCRR)
*
            nSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
*
            ij = (mdcj-1)*ndc + mdci
*                                                                      *
************************************************************************
*                                                                      *
            Call Do_NInt(AOInt,nAOInt,mGrid,Scr,iCmp,iBas_Eff,
     &                   TabAO_Pack(ipTabAO(jList_s):),jCmp,jBas_Eff,
     &                   nGrid_Tot,iSpin,mAO,nFn)
*
            Do iD = 1, iSpin
*
               If (nIrrep.ne.1) Then
*
*                 Gather the proper set of integral in the correct order
*                 for the symmetry adaptation.
*
*---------------- Distribute contributions to the SO integrals
*
                  Call SymAdd2(iAng,jAng,iCmp,jCmp,
     &                         iShell,jShell,iShll,jShll,
     &                         iAO,jAO,AOInt(1,iD),
     &                         iBas,iBas_Eff,jBas,jBas_Eff,
     &                         nOp,
     &                         iSkal,jSkal,Fact(ij),
     &                         FckInt(1,iD),nFclInt)
*
               Else
*
*                 Here we simply accumulate the AO integrals to the
*                 correct position in the full FckInt matrix.
*                 This should eventually be done outside of the loop.
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
         End Do                     ! jlist_s
      End Do                        ! ilist_s
*                                                                      *
************************************************************************
*                                                                      *
      Flop=Flop+DBLE(nGrid_Tot)
*
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(nSym)
#endif
      End
#endif
