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
      SubRoutine M2Grd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,la,lb,A,RB,nHer,
     &                 Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                 IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                 iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of M2 integrals used in   *
*         ECP calculations. The operator is a s-type gaussian          *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              DCR                                                     *
*              CrtCmp                                                  *
*              Assmbl                                                  *
*              CmbnMP                                                  *
*              DaXpY   (ESSL)                                          *
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
************************************************************************
      use Basis_Info
      use Center_Info
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), TC(3), C(3),
     &       Array(nZeta*nArr), Ccoor(3), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
      Integer iStabM(0:nStabM-1), iDCRT(0:7), lOper(nComp),
     &          iuvwx(4), kOp(2), lOp(4),
     &          IndGrd(3,2), JndGrd(3,4)
      Logical ABeq(3), IfGrad(3,2), JfGrad(3,4), TstFnc, TF, EQ
*
*-----Statement function for Cartesian index
*
      nElem(k)=(k+1)*(k+2)/2
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
*
      iRout = 122
      iPrint = nPrint(iRout)
*     Call QEnter('M2Grd')
*
      iIrrep = 0
      iuvwx(1) = nStab(mdc)
      iuvwx(2) = nStab(ndc)
      lOp(1) = kOp(1)
      lOp(2) = kOp(2)
      nDAO = nElem(la)*nElem(lb)
*
      nip = 1
      ipA = nip
      nip = nip + nZeta
      ipB = nip
      nip = nip + nZeta
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+2)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+2)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer
      ipQxyz = nip
      nip = nip + nZeta*3*(la+2)*(lb+2)
      ipK = nip
      nip = nip + nZeta
      ipZ = nip
      nip = nip + nZeta
      ipPx= nip
      nip = nip + nZeta
      ipPy= nip
      nip = nip + nZeta
      ipPz= nip
      nip = nip + nZeta
      If (nip-1.gt.nArr*nZeta) Then
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Call ErrTra
         Write (6,*) ' Abend in M2Grd'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In M2Grd: A',' ',A,1,3)
         Call RecPrt(' In M2Grd: RB',' ',RB,1,3)
         Call RecPrt(' In M2Grd: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In M2Grd: Kappa',' ',rKappa,nAlpha,nBeta)
         Call RecPrt(' In M2Grd: Zeta',' ',Zeta,nAlpha,nBeta)
         Call RecPrt(' In M2Grd: P',' ',P,nZeta,3)
         Write (6,*) ' In M2Grd: la,lb,nHer=',la,lb,nHer
      End If
*
      iStrt = ipA
      Do 20 iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(iStrt),1)
         iStrt = iStrt + nAlpha
 20   Continue
*
      iStrt = ipB
      Do 21 iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(iStrt),nAlpha)
         iStrt = iStrt + 1
 21   Continue
*
*-----Loop over nuclear centers
*
      kdc=0
      Do 100 kCnttp = 1, nCnttp
         If (.Not.dbsc(kCnttp)%ECP) Go To 111
         If (dbsc(kCnttp)%nM2.eq.0) Go To 111
*
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab, nStab(kdc+kCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
            iuvwx(3) = nStab(kdc+kCnt)
            iuvwx(4) = nStab(kdc+kCnt)
*
            Do 102 lDCRT = 0, nDCRT-1
               lOp(3) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
               lOp(4) = lOp(3)
               Call OA(iDCRT(lDCRT),C,TC)
               If (EQ(A,RB).and.EQ(A,TC)) Go To 102
*
               Do 1011 iM2xp = 1, dbsc(kCnttp)%nM2
                  Gamma = dbsc(kCnttp)%M2xp(iM2xp)
                  If (iPrint.ge.99) Write (6,*) ' Gamma=',Gamma
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
*
                  nDisp = IndDsp(kdc+kCnt,iIrrep)
                  Do 220 iCar = 0, 2
                     JfGrad(iCar+1,3) = .False.
                     iCmp = 2**iCar
                     If ( TF(kdc+kCnt,iIrrep,iCmp) .and.
     &                    .Not.dbsc(kCnttp)%pChrg ) Then
                        nDisp = nDisp + 1
                        If (Direct(nDisp)) Then
*--------------------------Reset flags for the basis set centers so that
*                          we will explicitly compute the derivatives
*                          with respect to those centers. Activate flag
*                          for the third center so that its derivative
*                          will be computed by the translational
*                          invariance.
                           JndGrd(iCar+1,1) = Abs(JndGrd(iCar+1,1))
                           JndGrd(iCar+1,2) = Abs(JndGrd(iCar+1,2))
                           JndGrd(iCar+1,3) = -nDisp
                           JfGrad(iCar+1,1) = .True.
                           JfGrad(iCar+1,2) = .True.
                        Else
                           JndGrd(iCar+1,3) = 0
                        End If
                     Else
                        JndGrd(iCar+1,3) = 0
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
*-----------------Modify the original basis.
*
                  Do 1012 iZeta = 1, nZeta
                     PTC2 = (P(iZeta,1)-TC(1))**2
     &                    + (P(iZeta,2)-TC(2))**2
     &                    + (P(iZeta,3)-TC(3))**2
                     Tmp0 = Zeta(iZeta)+Gamma
                     Tmp1 = Exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
                     Array(ipK+iZeta-1)    = rKappa(iZeta) * Tmp1
                     Array(ipZ+iZeta-1)    = Tmp0
                     Array(ipPx+iZeta-1)   =
     &                  (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
                     Array(ipPy+iZeta-1)   =
     &                  (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
                     Array(ipPz+iZeta-1)   =
     &                  (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
 1012             Continue
                  If (iPrint.ge.99) Then
                     Write (6,*) ' The modified basis set'
                     Call RecPrt(' In M2Grd: Kappa',' ',
     &                            Array(ipK),nAlpha,nBeta)
                     Call RecPrt(' In M2Grd: Zeta',' ',
     &                            Array(ipZ),nAlpha,nBeta)
                     Call RecPrt(' In M2Grd: P',' ',Array(ipPx),nZeta,3)
                     Call RecPrt(' In M2Grd: TC',' ',TC,1,3)
                  End If
*
*-----------------Compute the cartesian values of the basis functions
*                 angular part
*
                  ABeq(1) = A(1).eq.RB(1) .and. A(1).eq.TC(1)
                  ABeq(2) = A(2).eq.RB(2) .and. A(2).eq.TC(2)
                  ABeq(3) = A(3).eq.RB(3) .and. A(3).eq.TC(3)
                  Call CrtCmp(Array(ipZ),Array(ipPx),nZeta,A,
     &                        Array(ipAxyz),la+1,HerR(iHerR(nHer)),
     &                        nHer,ABeq)
                  Call CrtCmp(Array(ipZ),Array(ipPx),nZeta,RB,
     &                        Array(ipBxyz),lb+1,HerR(iHerR(nHer)),
     &                        nHer,ABeq)
*
*-----------------Compute the contribution from the multipole moment
*                 operator
*
                  ABeq(1) = .False.
                  ABeq(2) = .False.
                  ABeq(3) = .False.
                  Call CrtCmp(Array(ipZ),Array(ipPx),nZeta,Ccoor,
     &                        Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &                        nHer,ABeq)
*
*-----------------Compute the cartesian components for the multipole
*                 moment integrals. The integrals are factorized into
*                 components.
*
                  Call Assmbl(Array(ipQxyz),
     &                        Array(ipAxyz),la+1,
     &                        Array(ipRxyz),nOrdOp,
     &                        Array(ipBxyz),lb+1,
     &                        nZeta,HerW(iHerW(nHer)),nHer)
*
*-----------------Combine the cartesian components to the full one
*                 electron integral gradient.
*
                  Factor = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M2cf(iM2xp)
     &                   * Fact
                  Call CmbnM2(Array(ipQxyz),nZeta,la,lb,
     &                        Array(ipZ),Array(ipK),Final,
     &                        Array(ipA),Array(ipB),JfGrad,Factor,mVec)
                  If (iPrint.ge.99) Call RecPrt(' Final in M2Grd',' ',
     &                Final,nZeta*nElem(la)*nElem(lb),mVec)
*
*-----------------Distribute the gradient contributions
*
                  Call DistG1X(Final,DAO,nZeta,nDAO,mVec,Grad,nGrad,
     &                         JfGrad,JndGrd,iuvwx,lOp,iChBas,MxFnc,
     &                         nIrrep)
*
 1011          Continue
*
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
*
 100  Continue
*
*
*     Call QExit('M2Grd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(lOper)
      End If
      End
