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
* Copyright (C) 1991,2021, Roland Lindh                                *
************************************************************************
      Subroutine SymAdp_Full(lOper,
     &                       AOIntegrals, nBfn, PrpInt, nPrp,
     &                       list_s,nlist_s,Fact,ndc)
************************************************************************
*                                                                      *
* Object: to transform the one-electon matrix elements from AO basis   *
*         to SO basis.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1991                                             *
************************************************************************
      use iSD_data
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info,     only: iAOtSO
      use nq_Grid,       only: iBfn_Index
      use Basis_Info,    only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOIntegrals(nBfn,nBfn), PrpInt(nPrp), Fact(ndc,ndc)
      Integer list_s(2,nlist_s)
      Integer nOp(2)
      Integer, Parameter:: iTwoj(0:7)=[1,2,4,8,16,32,64,128]
      Integer :: jIC(0:7)=[-99,-99,-99,-99,-99,-99,-99,-99]
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      iIC = 1
      Do iIrrep = 0, nIrrep-1
         If (iAnd(lOper,iTwoj(iIrrep)).eq.0) Cycle
         jIC(iIrrep) = iIC
         iIC = iIC + 1
      End Do
*
      loper=1
      Do j1 = 0, nIrrep-1
         xa = DBLE(iChTbl(j1,nOp(1)))
         Do iBfn = 1, nBfn
            ilist_s = iBfn_Index(2,iBfn)
            iCmp    = iBfn_Index(3,iBfn)
            iB_Eff  = iBfn_Index(4,iBfn)
            indAO1  = iBfn_Index(6,iBfn)
            iSkal   = list_s(1,ilist_s)
            kDCRE   = list_s(2,ilist_s)
            iAO     = iSD( 7,iSkal)
            mdci    = iSD(10,iSkal)
            iShell  = iSD(11,iSkal)
            nOp(1) = NrOpr(kDCRE)
            If (iAOtSO(iAO+i1,j1)<0) Cycle
            iSO1=iAOtSO(iAO+i1,j1)

            Do j2 = 0, nIrrep-1
               j12 = iEor(j1,j2)
               If (iAnd(lOper,iTwoj(j12)).eq.0) Cycle

               iPnt = iPntSO(j1,j2,lOper,nbas)
               xb = DBLE(iChTbl(j2,nOp(2)))

               Do jBfn = 1, nBfn
                  jlist_s = iBfn_Index(2,jBfn)
                  jCmp    = iBfn_Index(3,jBfn)
                  iB_Eff  = iBfn_Index(4,iBfn)
                  indAO2  = iBfn_Index(6,iBfn)
                  jSkal   = list_s(1,jlist_s)
                  kDCRR   = list_s(2,jlist_s)
                  jAO     = iSD( 7,jSkal)
                  mdcj    = iSD(10,jSkal)
                  jShell  = iSD(11,jSkal)
                  nOp(2) = NrOpr(kDCRR)
                  If (iAOtSO(jAO+i2,j2)<0) Cycle
                  iSO2=iAOtSO(jAO+i2,j2)
*
                  iSO=iSO1+IndAO1-1
                  jSO=iSO2+IndAO2-1

*                 Diagonal block. Store only unique elements
                  If (.NOT.(j1.eq.j2 .and. iSO1.eq.iSO2 .and.
     &                      iSO<jSO)) Then

                     If (j1.eq.j2) Then
*------------           Diagonal symmetry block
                        Indij=iPnt + iTri(iSO,jSO)
                     Else
*------------           Off-diagonal symmetry block j1>j2
                        nRow = nBas(j1)
                        Indij=iPnt + nRow*(jSO-1)*nRow + iSO
                     End If

                     PrpInt(Indij) = PrpInt(Indij)
     &                             +Fact(mdci,mdcj)*xa*xb
     &                             *AOIntegrals(iBfn,jBfn)
                  End If
*
               End Do ! jBfn
            End Do    ! j2

         End Do       ! iBfn
      End Do          ! j1
*
      Return
      End
