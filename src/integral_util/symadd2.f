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
      Subroutine SymAdd2(iAng,jAng,iCmp,jCmp,iShell,jShell,
     &                   iShll,jShll,iAO,jAO,
     &                               AOInt,iBas,iBas_Eff,
     &                                     jBas,jBas_Eff,
     &                   nOp,iSkal,jSkal,Fact,PrpInt,nPrp)
************************************************************************
*                                                                      *
* Object: to transform the one-electon matrix elements from AO basis   *
*         to SO basis.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1991                                             *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 AOInt(iBas_Eff*jBas_Eff,iCmp,jCmp)
      Real*8 PrpInt(nPrp)
      Integer nOp(2)
      Integer iTwoj(0:7), jIC(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *

*
      iRout = 133
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Write (6,*) ' lOper=',lOper
         Call RecPrt(' In SymAdd: AOInt',' ',AOInt,iBas*jBas,
     &                iCmp*jCmp)
      End If

*... to be simplified

      loper=1
      jIC(:)=-999999999
      jIC(0) = 1
*
*     We denote the basis functions as X(alpha,mu,m_l,A)
*
*     alpha : irrep index
*     mu    : radial contraction index
*     m_l   : angular index
*     A     : Atomic center

      iAdd = iBas-iBas_Eff
      jAdd = jBas-jBas_Eff

      Do 100 j1 = 0, nIrrep-1
         xa = DBLE(iChTbl(j1,nOp(1)))
         Do 200 i1 = 1, iCmp
            If (iAOtSO(iAO+i1,j1)<0) Cycle
*
*           Do 300 j2 = 0, nIrrep-1
               j2=j1
               j12 = iEor(j1,j2)
*
*              If (iAnd(lOper,iTwoj(j12)).eq.0) Go To 300
               xb = DBLE(iChTbl(j2,nOp(2)))
*
*              If two sets of basis functions have the same
*              irrep index, and share atomic center and radial
*              contractions retrict the angular index combinations
*              to be the unique combinations.

*
               Do 400 i2 = 1, jCmp
                  If (iAOtSO(jAO+i2,j2)<0) Cycle

                  If (iShell.eq.jShell .and. j1.eq.j2 .and.
     &               i1<i2) Cycle

                  iSO1=iAOtSO(iAO+i1,j1)
                  iSO2=iAOtSO(jAO+i2,j2)

                  iPnt = iPntSO(j1,j2,lOper,nbas)
*
                  Do iB_Eff = 1, iBas_Eff
                     indAO1 = iB_Eff + iAdd
                     Do jB_Eff = 1, jBas_Eff
                        indAO2 = jB_Eff + jAdd

                        iSO=iSO1+IndAO1-1
                        jSO=iSO2+IndAO2-1
*
                        iFrom=(jB_Eff-1)*iBas_Eff+iB_Eff
*                       If (j1.eq.j2) Then
*------------            Diagonal symmetry block
                           If (iSO1.eq.iSO2 .and. iSO<jSO) Cycle
                           Indij=iPnt + iTri(iSO,jSO)
*                       Else
*------------              Off-diagonal symmetry block j1>j2
*                          nRow = nBas(j1)
*                          Indij=iPnt + nRow*(jSO-1)*nRow + iSO
*                       End If

                        PrpInt(Indij) = PrpInt(Indij)
     &                                +Fact*xa*xb*AOInt(iFrom,i1,i2)
                     End Do
                  End Do

                  If (.Not.(iShell.eq.jShell.and.nOp(1).ne.nOp(2)))
     &               Cycle

                  Do iB_Eff = 1, iBas_Eff
                     indAO1 = iB_Eff + iAdd
                     Do jB_Eff = 1, jBas_Eff
                        indAO2 = jB_Eff + jAdd

                        iSO=iSO1+IndAO2-1
                        jSO=iSO2+IndAO1-1

*           Diagonal block. Store only unique elements

*
                        iFrom=(jB_Eff-1)*iBas_Eff+iB_Eff
*                       If (j1.eq.j2) Then
*------------              Diagonal symmetry block
                           If (iSO1.eq.iSO2 .and. iSO<jSO) Cycle
                           Indij=iPnt + iTri(iSO,jSO)
*                       Else
*------------              Off-diagonal symmetry block j1>j2
*                          nRow = nBas(j1)
*                          Indij=iPnt + nRow*(jSO-1)*nRow + iSO
*                       End If

                        PrpInt(Indij) = PrpInt(Indij)
     &                                +Fact*xa*xb*AOInt(iFrom,i2,i1)
*
                     End Do
                  End Do
*
 400           Continue
*
*300        Continue
*
 200     Continue
 100  Continue

*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iAng)
         Call Unused_integer(jAng)
         Call Unused_integer(iShll)
         Call Unused_integer(jShll)
         Call Unused_integer(iSkal)
         Call Unused_integer(jSkal)
      End If
      End
