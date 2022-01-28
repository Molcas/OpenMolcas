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
* Copyright (C) 2000,2022, Roland Lindh                                *
*               Ajitha Devarajan                                       *
************************************************************************
      Subroutine DFT_Int(list_s,nlist_s,FckInt,nFckInt,nD,Fact,ndc,
     &                   list_bas)
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
      use nq_Grid, only: TabAO, Grid_AO, iBfn_Index,
     &                   AOIntegrals => Dens_AO
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
#include "stdalloc.fh"
      Real*8 Fact(ndc**2), FckInt(nFckInt,nD)
      Integer list_s(2,nlist_s), list_bas(2,nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the desired AO integrand here from the AOs, accumulate
*     contributions to the SO integrals on the fly.
*
      mGrid=SIZE(TabAO,2)
      mAO = SIZE(TabAO,1)
*
      nBfn = Size(AOIntegrals,1)
      nFn  = Size(Grid_AO,1)
      Call Do_NInt_d(mGrid,Grid_AO, TabAO,nBfn,nD,mAO,nFn)
      Call Do_NIntX(AOIntegrals,mGrid,Grid_AO,TabAO,nBfn,nD,mAO,nFn)
*                                                                      *
************************************************************************
*                                                                      *
*     Set up an indexation translation between the running index of
*     the AOIntegrals and the actual basis function index
*
      Call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
      iBfn_Index(:,:)=0

      iBfn = 0
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         iCmp  = iSD( 2,iSkal)
         nBas  = iSD( 3,iSkal)
         nBas_Eff=list_bas(1,ilist_s)
         iAO   = iSD( 7,iSkal)
         mdci  = iSD(10,iSkal)

         iAdd = nBas-nBas_Eff
         Do i1 = 1, iCmp
            iSO1 = iAOtSO(iAO+i1,0) ! just used when nIrrep=1
            Do i2 = 1, nBas_Eff
               IndAO1 = i2 + iAdd
               Indi = iSO1 + IndAO1 -1

               iBfn = iBfn + 1
               iBfn_Index(1,iBfn) = Indi
               iBfn_Index(2,iBfn) = ilist_s
               iBfn_Index(3,iBfn) = i1
               iBfn_Index(4,iBfn) = i2
               iBfn_Index(5,iBfn) = mdci
               iBfn_Index(6,iBfn) = IndAO1
            End Do
         End Do
      End Do
      If (iBfn.ne.nBfn) Then
         Write (6,*) 'Inconsistent numbers, iBfn/=nBfn:',iBfn,nBfn
         Call abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Do iD = 1, nD
         If (nIrrep.eq.1) Then
            Call AOAdd_Full(AOIntegrals(:,:,iD),nBfn,FckInt(:,iD),
     &                      nFckInt)
         Else
            Call SymAdp_Full(AOIntegrals(:,:,iD),nBfn,FckInt(:,iD),
     &                       nFckInt,list_s,nlist_s,Fact,ndc)
         End If
      End Do
      Call mma_deAllocate(iBfn_Index)
*                                                                      *
************************************************************************
*                                                                      *
      End
