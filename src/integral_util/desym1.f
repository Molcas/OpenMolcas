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
      Subroutine Desym1(lOper,iAng,jAng,iCmp,jCmp,iShell,jShell,
     &                  iShll,jShll,iAO,jAO,DAO,iBas,jBas,
     &                  DSO,nDSO,nOp,FactNd,Scrt)
************************************************************************
*                                                                      *
* Object: desymmetrize the first order density matrix.                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 DAO(iBas*jBas,iCmp,jCmp), DSO(iBas*jBas,nDSO),
     &       Scrt(iBas*jBas)
      Integer nOp(2)
*
      iRout = 133
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Write (6,*) ' lOper=',lOper
         Call RecPrt(' In Desym1: DSO',' ',DSO,iBas*jBas,nDSO)
      End If
*
      call dcopy_(iBas*jBas*iCmp*jCmp,[Zero],0,DAO,1)
*
*     D(P,Q)_ij = Sum(iSym,jSym) X(iSym,P) X(jSym,Q) D(iSym,jSym)_ij
*
*     Loop over DSO, iSym>=jSym
*
      lSO = 0
      Do 100 j1 = 0, nIrrep-1
         Xa= DBLE(iChTbl(j1,nOp(1)))
         Do 200 i1 = 1, iCmp
            If (iAOtSO(iAO+i1,j1)<0) Cycle
*
            Do 300 j2 = 0, j1
               j12 = iEor(j1,j2)
               If (iAnd(lOper,2**j12).eq.0) Go To 300
               Xb = DBLE(iChTbl(j2,nOp(2)))
               jMx = jCmp
               If (iShell.eq.jShell .and. j1.eq.j2) jMx = i1
               Do 400 i2 = 1, jMx
                  If (iAOtSO(jAO+i2,j2)<0) Cycle
                  lSO = lSO + 1
*
                  Deg=Two
                  If (j1.eq.j2) Deg=One
*
*-----------------Parity factor due to symmetry operations applied to
*                 angular part of the basis function.
*
                  Call DaXpY_(iBas*jBas,Deg*Xa*Xb,
     &                       DSO(1,lSO),1,
     &                       DAO(1,i1,i2),1)
*
                  If (iShell.eq.jShell .and. j1.eq.j2 .and.
     &                i1.ne.i2) Then
                     Call DGeTMO(DSO(1,lSO),iBas,iBas,jBas,Scrt,jBas)
                     Call DaXpY_(iBas*jBas,Deg*Xa*Xb,
     &                          Scrt,1,
     &                          DAO(1,i2,i1),1)
                  End If
 400           Continue
 300        Continue
*
 200     Continue
 100  Continue
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In Desym1: DAO',' ',DAO,iBas*jBas,iCmp*jCmp)
      End If
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iAng)
         Call Unused_integer(jAng)
         Call Unused_integer(iShll)
         Call Unused_integer(jShll)
         Call Unused_real(FactNd)
      End If
      End
