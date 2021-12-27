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
     &                       AOIntegrals, nBfn, Prpnt, nPrp,
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
      Real*8 AOIntIntegral(nBfn,nBfn), PrpInt(nPrp), Fact(ndc**2)
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
            i1      = iBfn_Index(3,iBfn)
            indAO1  = iBfn_Index(4,iBfn)
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
               kIC = jIC(j12)
               xb = DBLE(iChTbl(j2,nOp(2)))

               Do jBfn = 1, nBfn
                  jlist_s = iBfn_Index(2,jBfn)
                  i2      = iBfn_Index(3,jBfn)
                  indAO1  = iBfn_Index(4,iBfn)
                  jSkal   = list_s(1,jlist_s)
                  kDCRR   = list_s(2,jlist_s)
                  jAO     = iSD( 7,jSkal)
                  mdcj    = iSD(10,jSkal)
                  jShell  = iSD(11,jSkal)
                  nOp(2) = NrOpr(kDCRR)
                  If (iAOtSO(jAO+i2,j2)<0) Cycle
                  iSO2=iAOtSO(jAO+i2,j2)
*
                  ij = (mdcj-1)*ndc + mdci

*                 For the case of two identical shells and
*                 a diagonal irrep block restrict the loop
*                 to be triangular in the block. We achieve
*                 this by looping over the unique angular
*                 combinations. For the diagonal term of the
*                 angular indicies we will allow a redundance
*                 in the radial index. This was later compensated
*                 in SOAdd, which has a strict triangularization
*                 here too.
*
                  If (iShell.eq.jShell .and. j1.eq.j2
     &                .and. i2>i1) Cycle

                  If (j1.eq.j2.and.iSO1.lt.iSO2) Cycle

                  If (j1.eq.j2) Then
*------------        Diagonal symmetry block

                     Indij=iPnt + iTri(iSO1+indAO1,iSO2+indAO2)
                  Else
*------------        Off-diagonal symmetry block j1>j2

                     Indi = iSO1+indAO1
                     Indj = iSO2+indAO2
                     nRow = nBas(j1)
                     Indij = iPnt + nRow*(Indj-1) + Indi
                  End If

                  PrpInt(Indij) = PrpInt(Indij)
     &                          + Factor(ij)*xa*xb
     &                          * AOIntegral(iBfn,jBfn)

*                 If (.NOT.(iShell.eq.jShell.and.nOp(1).ne.nOp(2)))Cycle

*                 PrtInt(b,a) = PrtInt(b,a)
*    &                        + Factor(ij)*xa*xb
*    &                        * AOIntegral(iBfn,jBfn)
*                 End If

               End Do
            End Do

         End Do
      End Do
*
      Return
      End
