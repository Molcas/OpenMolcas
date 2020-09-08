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
      SubRoutine DrvN1_mck(Grad,nGrad)
************************************************************************
*                                                                      *
* Object: to compute the molecular gradient contribution due to the    *
*         nuclear repulsion energy.                                    *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October 1991                                             *
************************************************************************
      use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "disp2.fh"
      Real*8 A(3), B(3), RB(3), Grad(nGrad)
      Integer iDCRR(0:7)
      Logical EQ, TstFnc
*
      iRout = 33
      iPrint = nPrint(iRout)
      Call qEnter('DrvN1')
*
      iIrrep = 0
      mdc = 0
*-----Loop over centers with the same change
      Do 100 iCnttp = 1, nCnttp
         ZA = dbsc(iCnttp)%Charge
         If (ZA.eq.Zero) Go To 101
*--------Loop over all unique centers of this group
         Do 110 iCnt = 1, dbsc(iCnttp)%nCntr
            A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
            ndc = 0
            Do 200 jCnttp = 1, iCnttp
               ZAZB = ZA * dbsc(jCnttp)%Charge
               If (ZAZB.eq.Zero) Go To 201
               jCntMx = dbsc(jCnttp)%nCntr
               If (iCnttp.eq.jCnttp) jCntMx = iCnt
               Do 210 jCnt = 1, jCntMx
                  B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
                  Fact = One
*                 Factor due to resticted summation
                  If (EQ(A,B)) Fact = Half
*
*                 Find the DCR for the two centers
*
                  Call DCR(LmbdR,iOper,nIrrep,
     &                     dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,
     &                     dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,
     &                     iDCRR,nDCRR)
*
                  PreFct = Fact*ZAZB*DBLE(nIrrep)/DBLE(LmbdR)
                  Do 300 iR = 0, nDCRR-1
                     Call OA(iDCRR(iR),B,RB)
                     nOp = NrOpr(iDCRR(iR),iOper,nIrrep)
                     If (EQ(A,RB)) Go To 300
                     r12 = Sqrt((A(1)-RB(1))**2 +
     &                          (A(2)-RB(2))**2 +
     &                          (A(3)-RB(3))**2 )
*
*                    The factor u/g will ensure that the value of the
*                    gradient in symmetry adapted and no symmetry basis
*                    will have the same value.
*
                     nDisp = IndDsp(mdc+iCnt,iIrrep)
                     igu=nIrrep/dc(mdc+iCnt)%nStab
                     Do 400 iCar = 0, 2
                        iComp = 2**iCar
                        If ( TstFnc(iOper,nIrrep,iCoSet(0,0,mdc+iCnt),
     &                     iChTbl,iIrrep,iComp,dc(mdc+iCnt)%nStab)
     &                     ) Then
                           nDisp = nDisp + 1
                           If (.true.) Grad(nDisp) =
     &                        Grad(nDisp) - One/DBLE(igu) *
     &                        PreFct*(A(iCar+1)-RB(iCar+1))/(r12**3)
                        End If
 400                 Continue
*
                     nDisp = IndDsp(ndc+jCnt,iIrrep)
                     igv=nIrrep/dc(ndc+jCnt)%nStab
                     Do 450 iCar = 0, 2
                        iComp = 2**iCar
                        If ( TstFnc(iOper,nIrrep,iCoSet(0,0,ndc+jCnt),
     &                     iChTbl,iIrrep,iComp,dc(ndc+jCnt)%nStab)
     &                     ) Then
                           nDisp = nDisp + 1
                           If (.true.) Then
                              ps = DBLE(iPrmt(nOp,iChBas(2+iCar)))
                              Grad(nDisp) = Grad(nDisp) + ps *
     &                           One/DBLE(igv) *
     &                           PreFct*(A(iCar+1)-RB(iCar+1))/(r12**3)
                           End If
                        End If
 450                 Continue
 300              Continue
*
 210           Continue
 201           Continue
               ndc = ndc + dbsc(jCnttp)%nCntr
 200        Continue
 110     Continue
 101     Continue
         mdc = mdc + dbsc(iCnttp)%nCntr
 100  Continue
*
      Call qExit('DrvN1')
      Return
      End
