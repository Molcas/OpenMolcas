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
      SubRoutine Distg1X(g1,PAO,nT,mPAO,mVec,Grad,nGrad,IfGrad,IndGrd,
     &                   iStab,kOp)
************************************************************************
*                                                                      *
* Object: trace the gradient of the ERI's with the second order        *
*         density matrix                                               *
*                                                                      *
* Called from: Rysg1                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DGeMV   (ESSL)                                          *
*              DCopy   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChBas
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 g1(nT,mPAO,mVec), PAO(nT,mPAO), Grad(nGrad),
     &       Temp(9), PAOg1(12), Prmt(0:7)
      Logical IfGrad(3,4)
      Integer   IndGrd(3,4), kOp(4), iStab(4)
#ifdef _DEBUG_
      Character*80 Label
#endif
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
#ifdef _DEBUG_
      iRout = 239
      iPrint = nPrint(iRout)
      Call qEnter('Distg1')
      If (iPrint.ge.99) Then
         Call RecPrt('PAO',' ',PAO,nT,mPAO)
         Do 500 iVec = 1, mVec
            Write (Label,'(A,I2,A)') ' g1(',iVec,')'
            Call RecPrt(Label,' ',g1(1,1,iVec),nT,mPAO)
 500     Continue
         Call RecPrt('Accumulated gradient on entrance',
     &               ' ',Grad,nGrad,1)
      End If
      If (iPrint.ge.49) Write (6,*) IndGrd
#endif
*
*-----Trace the integral derivatives with the second order density
*     matrix.
*
      Call dGeMV_('T',nT*mPAO,mVec,
     &           One,g1,nT*mPAO,
     &           PAO,1,
     &           Zero,Temp,1)
      nVec = 0
#ifdef __INTEL_COMPILER
      Do kl = 1, 12
         iCar = (kl-1)/4 + 1
         iCent = kl - (iCar-1)*4
         ij = 3*(iCent-1)+iCar
         If (IfGrad(iCar,iCent)) Then
            nVec = nVec + 1
            PAOg1(ij) = Temp(nVec)
         Else
            PAOg1(ij) = Zero
         End If
      End Do
#else
*
*     Original code didn't work for Intel compiler with -O3
*     options since it swaps the loops.
*
      Do iCar = 1, 3
         Do iCent = 1, 4
            ij = 3*(iCent-1)+iCar
            If (IfGrad(iCar,iCent)) Then
               nVec = nVec + 1
               PAOg1(ij) = Temp(nVec)
            Else
               PAOg1(ij) = Zero
            End If
         End Do
      End Do
#endif
*
*-----Compute some of the contributions via the translational invariance
*
      Do 200 iCn = 1, 4
         Do 210 iCar = 1, 3
            If (IndGrd(iCar,iCn).lt.0) Then
               ij = 3*(iCn-1) + iCar
               Do 220 jCn = 1, 4
                  If (iCn.eq.jCn) Go To 220
                  If (IfGrad(iCar,jCn)) Then
                     kl = 3*(jCn-1) + iCar
                     PAOg1(ij)=PAOg1(ij)-PAOg1(kl)
                  End If
 220           Continue
            End If
 210     Continue
 200  Continue
#ifdef _DEBUG_
      If (iPrint.ge.49) Call RecPrt('PAOg1',' ',PAOg1,12,1)
#endif
*
*-----Distribute contribution to the gradient.
*
      Do 100 iCn = 1, 4
         Do 110 iCar = 1, 3
            ij = 3*(iCn-1) + iCar
            If (IndGrd(iCar,iCn).ne.0) Then
               iGrad = Abs(IndGrd(iCar,iCn))
*--------------Parity due to integration direction
               ps = xPrmt(kOp(iCn),iChBas(1+iCar))
               Fact = ps * DBLE(iStab(iCn)) / DBLE(nIrrep)
               Grad(iGrad) = Grad(iGrad) + Fact * PAOg1(ij)
            End If
 110    Continue
 100  Continue
#ifdef _DEBUG_
      If (iPrint.ge.49) Then
         Call RecPrt('Accumulated gradient on exit',
     &               ' ',Grad,nGrad,1)
      End If
*
      Call qExit('Distg1')
#endif
      Return
      End
