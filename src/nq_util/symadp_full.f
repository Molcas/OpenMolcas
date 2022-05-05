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
! Copyright (C) 1991,2021, Roland Lindh                                *
!***********************************************************************
      Subroutine SymAdp_Full(SOIntegrals,nSOInt,list_s,nlist_s,Fact,ndc,&
     &                       nD)
!***********************************************************************
!                                                                      *
! Object: to transform the one-electon matrix elements from AO basis   *
!         to SO basis.                                                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1991                                             *
!***********************************************************************
      use iSD_data
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info,     only: iAOtSO
      use nq_Grid,       only: iBfn_Index
      use nq_Grid, only: AOIntegrals => Dens_AO
      use Basis_Info,    only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 SOIntegrals(nSOInt,nD), Fact(ndc,ndc)
      Integer list_s(2,nlist_s)
      Integer nOp(2)
      Integer, Parameter:: iTwoj(0:7)=[1,2,4,8,16,32,64,128]
      Integer, Allocatable:: BasList(:,:)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Statement functions
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j
!                                                                      *
!***********************************************************************
!                                                                      *
      nBfn=SIZE(iBfn_Index,2)
      Call mma_Allocate(BasList,2,nBfn,Label='BasList')
      loper=1
      Do j1 = 0, nIrrep-1
         iPnt = iPntSO(j1,j1,lOper,nbas)

         ! Pick up only basis functions which contribute to (j1,j1)
         mBfn=0
         Do iBfn = 1, nBfn
            ilist_s = iBfn_Index(2,iBfn)
            iCmp    = iBfn_Index(3,iBfn)
            indAO1  = iBfn_Index(6,iBfn)
            iSkal   = list_s(1,ilist_s)
            iAO     = iSD( 7,iSkal)
            iSO1=iAOtSO(iAO+iCmp,j1)
            If (iSO1<0) Cycle
            mBfn=mBfn+1
            BasList(1,mBfn)=iBfn
            BasList(2,mBfn)=iSO1+IndAO1-1
         End Do

         Do iBfn_ = 1, mBfn
            iBfn=BasList(1,iBfn_)
            iSO=BasList(2,iBfn_)

            ilist_s = iBfn_Index(2,iBfn)
            iSkal   = list_s(1,ilist_s)
            kDCRE   = list_s(2,ilist_s)
            mdci    = iSD(10,iSkal)
            iShell  = iSD(11,iSkal)
            nOp(1) = NrOpr(kDCRE)
            xa = DBLE(iChTbl(j1,nOp(1)))

            Do jBfn_= 1, iBfn_
               jBfn=BasList(1,jBfn_)
               jSO=BasList(2,jBfn_)

               jlist_s = iBfn_Index(2,jBfn)
               jSkal   = list_s(1,jlist_s)
               kDCRR   = list_s(2,jlist_s)
               mdcj    = iSD(10,jSkal)
               jShell  = iSD(11,jSkal)
               nOp(2) = NrOpr(kDCRR)
               xb = DBLE(iChTbl(j1,nOp(2)))

               xaxb=xa*xb
               If (iShell==jShell .and. nOp(1)/=nOp(2)                  &
     &             .and. iSO==jSO) xaxb=xaxb*Two

               Indij = iPnt + iTri(iSO,jSO)

               SOIntegrals(Indij,:) = SOIntegrals(Indij,:)              &
     &                       + Fact(mdci,mdcj)*xaxb                     &
     &                       * AOIntegrals(iBfn,jBfn,:)

            End Do ! jBfn
         End Do    ! iBfn
      End Do       ! j1
      Call mma_deAllocate(BasList)
!
      Return
      End
