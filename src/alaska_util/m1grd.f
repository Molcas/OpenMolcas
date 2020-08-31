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
* Copyright (C) 1993, Roland Lindh                                     *
*               1993, Per Boussard                                     *
************************************************************************
      SubRoutine M1Grd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,la,lb,A,RB,nRys,
     &                 Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                 IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                 iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of the M1 integrals used  *
*         ECP calculations. The operator is the nuclear attraction     *
*         operator times a s-type gaussian function.                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              DCR                                                     *
*              Rysg1                                                   *
*              QExit                                                   *
*                                                                      *
*      Alpha : exponents of bra gaussians                              *
*      nAlpha: number of primitives (exponents) of bra gaussians       *
*      Beta  : as Alpha but for ket gaussians                          *
*      nBeta : as nAlpha but for the ket gaussians                     *
*      Zeta  : sum of exponents (nAlpha x nBeta)                       *
*      ZInv  : inverse of Zeta                                         *
*      rKappa: gaussian prefactor for the products of bra and ket      *
*              gaussians.                                              *
*      P     : center of new gaussian from the products of bra and ket *
*              gaussians.                                              *
*      Final : array for computed integrals                            *
*      nZeta : nAlpha x nBeta                                          *
*      nComp : number of components in the operator (e.g. dipolmoment  *
*              operator has three components)                          *
*      la    : total angular momentum of bra gaussian                  *
*      lb    : total angular momentum of ket gaussian                  *
*      A     : center of bra gaussian                                  *
*      B     : center of ket gaussian                                  *
*      nRys  : order of Rys- or Hermite-Gauss polynomial               *
*      Array : Auxiliary memory as requested by ECPMem                 *
*      nArr  : length of Array                                         *
*      Ccoor : coordinates of the operator, zero for symmetric oper.   *
*      NOrdOp: Order of the operator                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
*             Physics, University of Stockholm, Sweden, October '93.   *
*                                                                      *
*             Modified to gradients, December '93 (RL).                *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, Cff2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), C(3),
     &       Array(nZeta*nArr), Ccoor(3), TC(3), CoorAC(3,2),
     &       Coori(3,4), Coora(3,4), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
      Integer iStabM(0:nStabM-1), iDCRT(0:7), lOper(nComp),
     &        iAnga(4), iuvwx(4), kOp(2), lOp(4),
     &        IndGrd(3,2), JndGrd(3,4)
      Logical EQ, IfGrad(3,2), JfGrad(3,4), TstFnc, TF
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
*
*     Call qEnter('M1Grd')
      iRout = 193
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In M1Grd: A',' ',A,1,3)
         Call RecPrt(' In M1Grd: RB',' ',RB,1,3)
         Call RecPrt(' In M1Grd: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In M1Grd: P',' ',P,nZeta,3)
         Write (6,*) ' In M1Grd: la,lb=',' ',la,lb
      End If
*
*-----Allocate Scratch for primitives and work area for HRR
*
      ip = 1
      ipA = ip
      ip = ip + nZeta
      ipB = ip
      ip = ip + nZeta
      ipDAO = ip
      ip = ip + nZeta*nElem(la)*nElem(lb)
      ipK = ip
      ip = ip + nZeta
      ipZ = ip
      ip = ip + nZeta
      ipZI = ip
      ip = ip + nZeta
      ipPx = ip
      ip = ip + nZeta
      ipPy = ip
      ip = ip + nZeta
      ipPz = ip
      ip = ip + nZeta
      If (ip-1.gt.nArr*nZeta) Then
         Write (6,*) ' ip-1.gt.nArr*nZeta (M1 section)'
         Write (6,*) ' nArr,nZeta=',nArr,nZeta
         Call Abend()
      End If
      nArray = nArr*nZeta - ip + 1
*
      iIrrep = 0
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
      call dcopy_(3, A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(3, A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      If (la.ge.lb) Then
         call dcopy_(3, A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
      iuvwx(1) = nStab(mdc)
      iuvwx(2) = nStab(ndc)
      lOp(1) = kOp(1)
      lOp(2) = kOp(2)
*
      ipAOff = ipA
      Do 200 iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
         ipAOff = ipAOff + nAlpha
 200  Continue
*
      ipBOff = ipB
      Do 210 iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
         ipBOff = ipBOff + 1
 210  Continue
*
*-----Loop over nuclear centers.
*
      kdc = 0
      Do 100 kCnttp = 1, nCnttp
         If (.Not.dbsc(kCnttp)%ECP) Go To 111
         If (dbsc(kCnttp)%nM1.eq.0) Go To 111
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
            iuvwx(3) = nStab(kdc+kCnt)
            iuvwx(4) = nStab(kdc+kCnt)
*
            Do 102 lDCRT = 0, nDCRT-1
               lOp(3) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
               lOp(4) = lOp(3)
               TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
               TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
               TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
*--------------Branch out if one-center integral
               If (EQ(A,RB).and.EQ(A,TC)) Go To 102
               If (iPrint.ge.99) Call RecPrt(' In M1Grd: TC',' ',TC,1,3)
               call dcopy_(3,A,1,Coora(1,1),1)
               call dcopy_(3,RB,1,Coora(1,2),1)
               call dcopy_(6,Coora(1,1),1,Coori(1,1),1)
               If (.Not.EQ(A,RB) .or. .Not.EQ(A,TC)) Then
                  Coori(1,1) = Coori(1,1)+One
*                 Coora(1,1) = Coora(1,1)+One
               End If
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coori(1,3),1)
               call dcopy_(3,TC,1,Coori(1,4),1)
               call dcopy_(3,TC,1,Coora(1,3),1)
               call dcopy_(3,TC,1,Coora(1,4),1)
*
               Do 1011 iM1xp=1, dbsc(kCnttp)%nM1
                  Gamma = dbsc(kCnttp)%M1xp(iM1xp)
*
                  Call ICopy(6,IndGrd,1,JndGrd,1)
                  Do 10 i = 1, 3
                     Do 11 j = 1, 2
                        JfGrad(i,j) = IfGrad(i,j)
 11                  Continue
 10               Continue
*
*-----------------Derivatives with respect to the operator is computed
*                 via the translational invariance.
*                 Some extra care is needed here due to that Rys2Dg will
*                 try to avoid some of the work.
*
                  nDisp = IndDsp(kdc+kCnt,iIrrep)
                  Do 220 iCar = 0, 2
*--------------------No direct assembly of contribution from the operat.
                     JfGrad(iCar+1,3) = .False.
                     JndGrd(iCar+1,3) = 0
                     iCmp = 2**iCar
                     If ( TF(kdc+kCnt,iIrrep,iCmp) .and.
     &                    .Not.pChrg(kCnttp) ) Then
*-----------------------Displacement is symmetric
                        nDisp = nDisp + 1
                        If (Direct(nDisp)) Then
*--------------------------Reset flags for the basis set centers so that
*                          we will explicitly compute the derivatives
*                          with respect to those centers. Activate flag
*                          for the third center so that its derivative
*                          will be computed by the translational
*                          invariance.
                           JfGrad(iCar+1,1) = .True.
                           JfGrad(iCar+1,2) = .True.
                           If ( A(iCar+1).ne.TC(iCar+1) .and.
     &                         RB(iCar+1).ne.TC(iCar+1)) Then
*-----------------------------Three center case
                              JndGrd(iCar+1,1) = Abs(JndGrd(iCar+1,1))
                              JndGrd(iCar+1,2) = Abs(JndGrd(iCar+1,2))
                              JndGrd(iCar+1,3) = -nDisp
                              JfGrad(iCar+1,1) = .True.
                              JfGrad(iCar+1,2) = .True.
                            Else If ( A(iCar+1).eq.TC(iCar+1) .and.
     &                               RB(iCar+1).ne.TC(iCar+1)) Then
*-----------------------------Two center case
                              JndGrd(iCar+1,1) = -Abs(JndGrd(iCar+1,1))
                              JndGrd(iCar+1,2) = Abs(JndGrd(iCar+1,2))
                              JfGrad(iCar+1,1) = .False.
                              JfGrad(iCar+1,2) = .True.
                            Else If ( A(iCar+1).ne.TC(iCar+1) .and.
     &                               RB(iCar+1).eq.TC(iCar+1)) Then
*-----------------------------Two center case
                              JndGrd(iCar+1,1) = Abs(JndGrd(iCar+1,1))
                              JndGrd(iCar+1,2) = -Abs(JndGrd(iCar+1,2))
                              JfGrad(iCar+1,1) = .True.
                              JfGrad(iCar+1,2) = .False.
                            Else
*-----------------------------One center case
                              JndGrd(iCar+1,1) = 0
                              JndGrd(iCar+1,2) = 0
                              JfGrad(iCar+1,1) = .False.
                              JfGrad(iCar+1,2) = .False.
                            End If
                        End If
                     End If
 220              Continue
*-----------------No derivatives with respect to the fourth center.
                  Call ICopy(3,[0],0,JndGrd(1,4),1)
                  JfGrad(1,4) = .False.
                  JfGrad(2,4) = .False.
                  JfGrad(3,4) = .False.
                  mGrad = 0
                  Do 231 iCar = 1, 3
                     Do 232 i = 1, 2
                        If (JfGrad(iCar,i)) mGrad = mGrad + 1
 232                 Continue
 231              Continue
                  If (iPrint.ge.99) Write (6,*) ' mGrad=',mGrad
                  If (mGrad.eq.0) Go To 1011
*
*-----------------Modify the original basis. Observe that
*                 simplification due to A=B are not valid for the
*                 exponent index, eq. P-A=/=0.
*
                  Do 1012 iZeta = 1, nZeta
                     PTC2 = (P(iZeta,1)-TC(1))**2
     &                    + (P(iZeta,2)-TC(2))**2
     &                    + (P(iZeta,3)-TC(3))**2
                     Tmp0 = Zeta(iZeta)+Gamma
                     Tmp1 = Exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
                     Array(ipK +iZeta-1) = rKappa(iZeta) * Tmp1
                     Array(ipZ +iZeta-1) = Tmp0
                     Array(ipZI+iZeta-1) = One/Tmp0
                     Array(ipPx+iZeta-1) =
     &                  (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
                     Array(ipPy+iZeta-1) =
     &                  (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
                     Array(ipPz+iZeta-1) =
     &                  (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
 1012             Continue
*
*-----------------Modify the density matrix with the prefactor
*
                  Fact = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M1cf(iM1xp)*
     &                   (DBLE(nStabM) / DBLE(LmbdT)) * Two * Pi
                  nDAO = nElem(la)*nElem(lb)
                  Do 300 iDAO = 1, nDAO
                     Do 310 iZeta = 1, nZeta
                        Fac = Fact * Array(ipK+iZeta-1) *
     &                        Array(ipZI+iZeta-1)
                        ipDAOt = nZeta*(iDAO-1) + iZeta-1 + ipDAO
                        Array(ipDAOt)= Fac * DAO(iZeta,iDAO)
 310                 Continue
 300              Continue
                  If (iPrint.ge.99) Then
                     Write (6,*) ' Charge=',dbsc(kCnttp)%Charge
                     Write (6,*) ' Fact=',Fact
                     Write (6,*) ' IndGrd=',IndGrd
                     Write (6,*) ' JndGrd=',JndGrd
                     Call RecPrt('DAO*Fact',' ',Array(ipDAO),nZeta,nDAO)
                  End If
*
*-----------------Compute integrals with the Rys quadrature.
*
                  Call Rysg1(iAnga,nRys,nZeta,
     &                       Array(ipA),Array(ipB),[One],[One],
     &                       Array(ipZ),Array(ipZI),nZeta,[One],[One],1,
     &                       Array(ipPx),nZeta,TC,1,Coori,Coora,
     &                       CoorAC,Array(ip),nArray,
     &                       TNAI1,Fake,Cff2D,
     &                       Array(ipDAO),nDAO,Grad,nGrad,JfGrad,JndGrd,
     &                       lOp,iuvwx)
*
 1011          Continue
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
*
*
*     Call QExit('M1Grd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_real_array(Final)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
      End If
      End
