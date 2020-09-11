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
      Subroutine DesymD(lOper,iAng,jAng,iCmp,jCmp,iShell,jShell,
     &                  iShll,jShll,iAO,jAO,DAO,
     &                  iBas,jBas,DSO,nDSO,nOp,FactNd)
************************************************************************
*                                                                      *
* Object: desymmetrize the first order density matrix.                 *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy   (ESSL)                                          *
*              DYaX    (ESSL)                                          *
*              DaXpY   (ESSL)                                          *
*              DScal   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      Use Basis_Info
      use Symmetry_Info, only: iChTbl
      use SOAO_Info, only: iAOtSO
      use Real_Spherical, only: iSphCr
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
      Real*8 DAO(iBas*jBas,iCmp,jCmp), DSO(iBas*jBas,nDSO), Prmt(0:7)
      Integer nOp(2)
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iRout = 133
      iPrint = nPrint(iRout)
      iQ=0
      If (iPrint.ge.99) Then
         Write (6,*) ' lOper=',lOper
         Call RecPrt(' In DesymD: DSO',' ',DSO,iBas*jBas,nDSO)
      End If
*
      call dcopy_(iBas*jBas*iCmp*jCmp,[Zero],0,DAO,1)
      lSO = 1
      ii = iAng*(iAng+1)*(iAng+2)/6
      jj = jAng*(jAng+1)*(jAng+2)/6
      Do 100 j1 = 0, nIrrep-1
         Xa= DBLE(iChTbl(j1,nOp(1)))
         Do 200 i1 = 1, iCmp
            If (iAOtSO(iAO+i1,j1)<0) Cycle
            iChBs = iChBas(ii+i1)
            If (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
            pa = xPrmt(iOper(nOp(1)),iChBs)
*
            Do 300 j2 = 0, j1
               j12 = iEor(j1,j2)
               If (iAnd(lOper,2**j12).eq.0) Go To 300
               Xb = DBLE(iChTbl(j2,nOp(2)))
               jMx = jCmp
               If (iShell.eq.jShell .and. j1.eq.j2) jMx = i1
               Do 400 i2 = 1, jMx
                  If (iAOtSO(jAO+i2,j2)<0) Cycle
                  jChBs = iChBas(jj+i2)
                  If (Shells(jShll)%Transf)
     &                jChBs = iChBas(iSphCr(jj+i2))
*
                  Deg=Two
                  If (j1.eq.j2 .and. iShell.eq.jShell .and.
     &                i1.eq.i2) Deg=One
*
*-----------------Parity factor due to symmetry operations applied to
*                 angular part of the basis function.
                  FactNs = pa * xPrmt(iOper(nOp(2)),jChBs)
                  Call DaXpY_(iBas*jBas,Deg*Xa*Xb*FactNs,
     &                       DSO(1,lSO),1,
     &                       DAO(1,i1,i2),1)
                  lSO = lSO + 1
 400           Continue
 300        Continue
*
 200     Continue
 100  Continue
*
      If (FactNd.ne.One) Call DScal_(iBas*jBas*iCmp*jCmp,FactNd,DAO,1)
      If (iPrint.ge.99) Then
         Call RecPrt(' In DesymD: DAO',' ',DAO,iBas*jBas,iCmp*jCmp)
      End If
      Return
      End
