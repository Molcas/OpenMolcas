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
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "disp2.fh"
#include "WrkSpc.fh"
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
         If (Charge(iCnttp).eq.Zero) Go To 101
         ZA = Charge(iCnttp)
         ixyz = dbsc(iCnttp)%ipCntr
*--------Loop over all unique centers of this group
         Do 110 iCnt = 1, dbsc(iCnttp)%nCntr
            A(1) = Work(ixyz  )
            A(2) = Work(ixyz+1)
            A(3) = Work(ixyz+2)
*
            ndc = 0
            Do 200 jCnttp = 1, iCnttp
               If (Charge(jCnttp).eq.Zero) Go To 201
               ZAZB = ZA * Charge(jCnttp)
               jxyz = dbsc(jCnttp)%ipCntr
               jCntMx = dbsc(jCnttp)%nCntr
               If (iCnttp.eq.jCnttp) jCntMx = iCnt
               Do 210 jCnt = 1, jCntMx
                  B(1) = Work(jxyz  )
                  B(2) = Work(jxyz+1)
                  B(3) = Work(jxyz+2)
*
                  Fact = One
*                 Factor due to resticted summation
                  If (EQ(A,B)) Fact = Half
*
*                 Find the DCR for the two centers
*
                  Call DCR(LmbdR,iOper,nIrrep,
     &                     jStab(0,mdc+iCnt),nStab(mdc+iCnt),
     &                     jStab(0,ndc+jCnt),nStab(ndc+jCnt),
     &                     iDCRR,nDCRR)
*
                  PreFct = Fact*ZAZB*DBLE(nIrrep)/DBLE(LmbdR)
                  Do 300 iR = 0, nDCRR-1
                     RB(1) = DBLE(iPhase(1,iDCRR(iR)))*B(1)
                     RB(2) = DBLE(iPhase(2,iDCRR(iR)))*B(2)
                     RB(3) = DBLE(iPhase(3,iDCRR(iR)))*B(3)
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
                     igu=nIrrep/nStab(mdc+iCnt)
                     Do 400 iCar = 0, 2
                        iComp = 2**iCar
                        If ( TstFnc(iOper,nIrrep,
     &                     iCoSet(0,0,mdc+iCnt),
     &                     nIrrep/nStab(mdc+iCnt),iChTbl,iIrrep,
     &                     iComp,nStab(mdc+iCnt)) ) Then
                           nDisp = nDisp + 1
                           If (.true.) Grad(nDisp) =
     &                        Grad(nDisp) - One/DBLE(igu) *
     &                        PreFct*(A(iCar+1)-RB(iCar+1))/(r12**3)
                        End If
 400                 Continue
*
                     nDisp = IndDsp(ndc+jCnt,iIrrep)
                     igv=nIrrep/nStab(ndc+jCnt)
                     Do 450 iCar = 0, 2
                        iComp = 2**iCar
                        If ( TstFnc(iOper,nIrrep,
     &                     iCoSet(0,0,ndc+jCnt),
     &                     nIrrep/nStab(ndc+jCnt),iChTbl,iIrrep,
     &                     iComp,nStab(ndc+jCnt)) ) Then
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
                  jxyz = jxyz + 3
 210           Continue
 201           Continue
               ndc = ndc + dbsc(jCnttp)%nCntr
 200        Continue
            ixyz = ixyz + 3
 110     Continue
 101     Continue
         mdc = mdc + dbsc(iCnttp)%nCntr
 100  Continue
*
      Call qExit('DrvN1')
      Return
      End
