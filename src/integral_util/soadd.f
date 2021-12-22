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
      SubRoutine SOAdd(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,
     &                  iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
************************************************************************
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 SOInt(iBas*jBas,nSOInt), PrpInt(nPrp)
      Logical AeqB
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *

*
      iRout = 130
      iPrint = nPrint(iRout)
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt(' In SOAdd:SOInt',' ',SOInt,iBas*jBas,nSOInt)
      End If
#endif
*
      lSO = 0
      Do 100 j1 = 0, nIrrep-1
       Do 200 i1 = 1, iCmp
        If (iAOtSO(iAO+i1,j1)<0) Cycle
*
*       Scatter the SO's onto lower rectangular blocks and triangular
*       diagonal blocks.
*
        Do 300 j2 = 0, j1
         j12 = iEor(j1,j2)
         If (iAnd(lOper,2**j12).eq.0) Cycle
         jjMx = jCmp
         If (iShell.eq.jShell .and. j1.eq.j2) jjMx = i1
         Do 400 i2 = 1, jjMx
          If (iAOtSO(jAO+i2,j2)<0) Cycle
          lSO = lSO + 1
          iSO1=iAOtSO(iAO+i1,j1)
          iSO2=iAOtSO(jAO+i2,j2)
*
          iPnt = iPntSO(j1,j2,lOper,nbas)
          Do 500 indAO1 = 0, iBas-1
*         Diagonal block. Store only unique elements
           jBsMax = jBas - 1
           If (j1.eq.j2 .and. iSO1.eq.iSO2) jBsMax=indAO1
           Do 600 indAO2 = 0, jBsMax
            ip = indAO2*iBas + indAO1 + 1
*
*           Move one electron integral.
*
            If (j1.eq.j2) Then
*------------Diagonal symmetry block
             Indij=iPnt + iTri(iSO1+indAO1,iSO2+indAO2)
            Else
*------------Off-diagonal symmetry block j1>j2
             Indi = iSO1+indAO1
             Indj = iSO2+indAO2
             nRow = nBas(j1)
             Indij=iPnt + nRow*(Indj-1) + Indi
            End If
            PrpInt(Indij) = PrpInt(Indij) + SOInt(ip,lSO)
*
 600       Continue
 500      Continue
*
 400     Continue
 300    Continue
*
 200   Continue
 100  Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_logical(AeqB)
      End If
      End
