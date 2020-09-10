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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Subroutine SymAdp(iAng, iCmp, jCmp, kCmp, lCmp, Shijij,
     &                  iShll, iShell, IndShl, kOp, ijkl,
     &                  Aux,nAux,AOInt,SOInt,nSOInt,Done)
************************************************************************
*  Object: to transform the integrals in AO basis to symmetry adapted  *
*          orbitals , SO. This is done by accumulating the AO inte-    *
*          grals onto the SO integrals.                                *
*          Observe that one of the operator is the Unit operation      *
*          performed on center A. However, since we scramble the order *
*          we do not really know which center is which. However, the   *
*          Unit operator will always give a factor of one. Hence, this *
*          is a convenient way to do the symmetry transformation with- *
*          out having to know the order.                               *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*          This code is never executed in the no symmetry case!!!      *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              GetMem                                                  *
*              DnaXpY   (ESSL)                                         *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      use Basis_Info
      use Symmetry_Info, only: iChTbl
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Real*8 AOInt(ijkl,iCmp,jCmp,kCmp,lCmp),
     &       SOInt(ijkl,nSOInt), Aux(nAux)
      Logical Shij, Shkl, Shijij, Qij, Qkl, Qijij, Done
      Integer iAng(4), iShell(4), iShll(4), kOp(4), IndShl(4)
*     Local Array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer iTwoj(0:7)
      Real*8 Prmt(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      iOff(ixyz)  = ixyz*(ixyz+1)*(ixyz+2)/6
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iRout = 38
      iPrint = nPrint(iRout)
      Done=.False.
      k12=0
      k34=0
      ii = iOff(iAng(1))
      jj = iOff(iAng(2))
      kk = iOff(iAng(3))
      ll = iOff(iAng(4))
      MemSO2 = 1
      If (iPrint.ge.99) Then
         Call RecPrt(' In SymAdp: AOInt ',' ',
     &               AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
      End If
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      Do 100 i1 = 1, iCmp
         Do 101 j = 0, nIrrep-1
            iSym(j) = iAnd(IrrCmp(IndShl(1)+i1),iTwoj(j))
101      Continue
         jCmpMx = jCmp
         If (Shij) jCmpMx = i1
         iChBs = iChBas(ii+i1)
         If (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
         pEa = xPrmt(iOper(kOp(1)),iChBs)
         Do 200 i2 = 1, jCmpMx
            Do 201 j = 0, nIrrep-1
               jSym(j) = iAnd(IrrCmp(IndShl(2)+i2),iTwoj(j))
201         Continue
            jChBs = iChBas(jj+i2)
            If (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
            pRb = xPrmt(iOper(kOp(2)),jChBs) * pEa
            Qij = i1.eq.i2
            If (iShell(2).gt.iShell(1)) Then
               i12 = jCmp*(i1-1) + i2
            Else
               i12 = iCmp*(i2-1) + i1
            End If
            Do 300 i3 = 1, kCmp
               Do 301 j = 0, nIrrep-1
                  kSym(j) = iAnd(IrrCmp(IndShl(3)+i3),iTwoj(j))
301            Continue
               lCmpMx = lCmp
               If (Shkl) lCmpMx = i3
               kChBs = iChBas(kk+i3)
               If (Shells(iShll(3))%Transf)
     &            kChBs = iChBas(iSphCr(kk+i3))
               pTc = xPrmt(iOper(kOp(3)),kChBs) * pRb
               Do 400 i4 = 1, lCmpMx
                  Do 401 j = 0, nIrrep-1
                     lSym(j) = iAnd(IrrCmp(IndShl(4)+i4),iTwoj(j))
401               Continue
                  Qkl = i3.eq.i4
                  lChBs = iChBas(ll+i4)
                  If (Shells(iShll(4))%Transf)
     &               lChBs = iChBas(iSphCr(ll+i4))
                  pTSd= xPrmt(iOper(kOp(4)),lChBs) * pTc
                  If (iShell(4).gt.iShell(3)) Then
                     i34 = lCmp*(i3-1) + i4
                  Else
                     i34 = kCmp*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) Go To 400
                  Qijij = Shijij .and. i12.eq.i34
*
*      Loop over irreps which are spanned by the basis function.
*      Again, the loop structure is restricted to ensure unique
*      integrals.
*
       iAux = 0
       Do 110 j1 = 0, nIrrep-1
          If (iSym(j1).eq.0) Go To 110
          Xa = DBLE(iChTbl(j1,kOp(1))) * pTSd
          j2Max = nIrrep-1
          If (Shij .and. Qij) j2Max = j1
          Do 210 j2 = 0, j2Max
             If (jSym(j2).eq.0) Go To 210
             Xb = DBLE(iChTbl(j2,kOp(2))) * Xa
             j12 = iEor(j1,j2)
             If (Qijij) Then
                If (Shij .and. Qij) Then
                   k12 = j1*(j1+1)/2 + j2+1
                Else If (Shij) Then
                   k12 = nIrrep*j1 + j2+1
                Else If (iShell(1).gt.iShell(2)) Then
                   k12 = nIrrep*j1 + j2+1
                Else
                   k12 = nIrrep*j2 + j1+1
                End If
             End If
             Do 310 j3 = 0, nIrrep-1
                If (kSym(j3).eq.0) Go To 310
                j4 = iEor(j12,j3)
                If (lSym(j4).eq.0) Go To 310
                If (Shkl .and. Qkl .and. j4.gt.j3) Go To 310
                If (Qijij) Then
                   If (Shkl .and. Qkl) Then
                      k34 = j3*(j3+1)/2 + j4+1
                   Else If (Shkl) Then
                      k34 = nIrrep*j3 + j4+1
                   Else If (iShell(3).gt.iShell(4)) Then
                      k34 = nIrrep*j3 + j4+1
                   Else
                      k34 = nIrrep*j4 + j3+1
                   End If
                   If (Qijij .and. k34.gt.k12) Go To 310
                End If
                Xg = DBLE(iChTbl(j3,kOp(3))) * Xb
                iAux = iAux + 1
                Aux(iAux) = DBLE(iChTbl(j4,kOp(4))) * Xg
*
 310         Continue
 210      Continue
 110   Continue
*
       If (iPrint.ge.99) Call RecPrt(' Aux',' ',Aux,iAux,1)
       If (iAux.ne.0) Then
          If (iAux.ne.1) Then
             Call DNaXpY(iAux,ijkl,Aux,1,AOInt(1,i1,i2,i3,i4),1,0,
     &                   SOInt(1,MemSO2),1,ijkl)
          Else
             Call DaXpY_(ijkl,Aux(1),AOInt(1,i1,i2,i3,i4),1,
     &                  SOInt(1,MemSO2),1)
          End If
          MemSO2 = MemSO2 + iAux
       End If
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
*
*     Call RecPrt(' On exit from SymAdp: SOInt ',' ',SOInt,ijkl,nSOInt)
      Return
      End
