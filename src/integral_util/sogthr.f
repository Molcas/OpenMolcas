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
      SubRoutine SOGthr(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,
     &                  iCmp,jCmp,iShell,jShell,IndShl,JndShl,
     &                  AeqB,iAO,jAO)
************************************************************************
*                                                                      *
* Object: to gather elements, from the Fock or 1st order density matrix*
*         in SO-basis, which are associated with a shell pair.         *
*         OBSERVE that the matrix is folded to triangular form, i.e.   *
*         the off diagonal elements has twice there value. This in     *
*         order to reduce the summation to i>=j. However, this will    *
*         not be used for the diagonal blocks. Hence, these elements   *
*         will only be assigned half the value from the original       *
*         matrix.                                                      *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
      Real*8 SOInt(iBas*jBas,nSOInt), PrpInt(nPrp)
      Logical AeqB
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
      iRout = 130
      iPrint = nPrint(iRout)
*     Call qEnter('SOGthr')
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In SOGthr: PrpInt',' ',PrpInt,1,nPrp)
         Write (6,*) ' iAO, jAO=',iAO, jAO
         Write (6,*) ' iShell, jShell=',iShell, jShell
      End If
      lSO = 0
      Do 100 j1 = 0, nIrrep-1
       Do 200 i1 = 1, iCmp
        If (iAnd(IrrCmp(IndShl+i1),2**j1).eq.0) Go To 200
*
*       Gather the SO's from lower rectangular blocks and triangular
*       diagonal blocks.
*
        Do 300 j2 = 0, j1
         j12 = iEor(j1,j2)
         If (iAnd(lOper,2**j12).eq.0) Go To 300
         jjMx = jCmp
         If (iShell.eq.jShell .and. j1.eq.j2) jjMx = i1
         Do 400 i2 = 1, jjMx
          If (iAnd(IrrCmp(JndShl+i2),2**j2).eq.0) Go To 400
          lSO = lSO + 1
          iSO1=iAOtSO(iAO+i1,j1)
          iSO2=iAOtSO(jAO+i2,j2)
*
          iPnt = iPntSO(j1,j2,lOper,nbas)
*         Write (*,*)  iSO1, iSO2, iPnt, i1, j1, i2, j2
          Do 500 indAO1 = 0, iBas-1
*         Diagonal block (only unique elements).
           Do 600 indAO2 = 0, jBas - 1
            ipij = indAO2*iBas + indAO1 + 1
            Indi = iSO1+indAO1
            Indj = iSO2+indAO2
*
*           Move one element and unfold.
*
            Fact=Half
            If (Indi.eq.Indj) Fact=One
*
            If (j1.eq.j2) Then
               SOInt(ipij,lSO)=Fact*PrpInt(iPnt + iTri(Indi,Indj))
            Else
               nRow = nBas(j1)
               SOInt(ipij,lSO)=Fact*PrpInt(iPnt + nRow*(Indj-1)+Indi)
            End If
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
      If (iPrint.ge.99) Then
         Call RecPrt(' In SOGthr: SOInt',' ',SOInt,iBas*jBas,nSOInt)
      End If
      If (iPrint.ge.99) Call GetMem(' Exit SOGthr','CHECK','REAL',
     &                              iDum,iDum)
*     Call qExit('SOGthr')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(AeqB)
      End
