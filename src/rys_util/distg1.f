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
      SubRoutine Distg1(Temp,mVec,Grad,nGrad,IfGrad,IndGrd,
     &                  iStab,kOp)
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
      Real*8 Grad(nGrad), Temp(9), PAOg1(12), Prmt(0:7)
      Logical IfGrad(3,4)
      Integer IndGrd(3,4), kOp(4), iStab(4)
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iRout = 239
      iPrint = nPrint(iRout)
*     Call qEnter('Distg1')
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt('Accumulated gradient on entrance',
     &               ' ',Grad,nGrad,1)
         Write (6,*) ' kOp=',kOp
         Write (6,*) ' iStab=',iStab
         Call RecPrt('Distg1: Temp',' ',Temp,9,1)
      End If
      If (iPrint.ge.49) Write (6,*) IndGrd
#endif
*
*----- Distribute Temp in PAOg1
*
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
*-----a) Compute some of the contributions via the translational
*        invariance
*-----b) Distribute contribution to the gradient.
*
      Do iCn = 1, 4
         Do iCar = 1, 3
            ij = 3*(iCn-1) + iCar
*
*           a)
*
            If (IndGrd(iCar,iCn).lt.0) Then
               Do jCn = 1, 4
                  If (iCn.ne.jCn.and.IfGrad(iCar,jCn)) Then
                     kl = 3*(jCn-1) + iCar
                     PAOg1(ij)=PAOg1(ij)-PAOg1(kl)
                  End If
               End Do
            End If
*
*           b)
*
            If (IndGrd(iCar,iCn).ne.0) Then
               iGrad = Abs(IndGrd(iCar,iCn))
*--------------Parity due to integration direction
               ps = xPrmt(kOp(iCn),iChBas(1+iCar))
               Fact = ps * DBLE(iStab(iCn)) / DBLE(nIrrep)
               Grad(iGrad) = Grad(iGrad) + Fact * PAOg1(ij)
            End If
*
         End Do
      End Do
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.49) Then
         Call RecPrt('PAOg1',' ',PAOg1,12,1)
         Call RecPrt('Accumulated gradient on exit',
     &               ' ',Grad,nGrad,1)
      End If
#endif
*
*     Call qExit('Distg1')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(mVec)
      End
