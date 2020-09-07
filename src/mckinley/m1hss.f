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
* Copyright (C) 1991, Anders Bernhardsson                              *
*               1991, Roland Lindh                                     *
************************************************************************
      SubRoutine M1Hss(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,la,lb,A,RB,nRys,
     &                 Array,nArr,Ccoor,nOrdOp,Hess,nHess,
     &                 IfHss,IndHss,ifgrd,IndGrd,DAO,mdc,ndc,nOp,
     &                 lOper,nComp,iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: to compute the gradient of the nuclear attraction integrals. *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy  (ESSL)                                           *
*              ICopy                                                   *
*              Rysg1                                                   *
*              QExit                                                   *
*                                                                      *
*             Anders Bernhardsson & Roland Lindh,                      *
*             Dept. of Theoretical Chemistry, University               *
*             of Lund, SWEDEN.                                         *
*             October 1991                                             *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
c#include "print.fh"
#include "disp.fh"
#include "disp2.fh"
      Integer IndGrd(0:2,0:1,0:(nIrrep-1)),
     &          IndHss(0:1,0:2,0:1,0:2,0:(nIrrep-1)),
     &          nOp(2), lOper(nComp), iStabM(0:nStabM-1),
     &          iDCRT(0:7)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), C(3), TC(3),
     &       Array(nArr), Ccoor(3), Hess(nHess),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
       Logical IfHss(0:1,0:2,0:1,0:2),IfGrd(0:2,0:1), TstFnc, TF,
     &         EQ,IfG(0:3),Tr(0:3)
*
*     Local arrrays
*
      Real*8 Coora(3,4), Coori(3,4), CoorAC(3,2)
      Integer iAnga(4), JndGrd(0:2,0:3,0:7),
     &           JndHss(0:3,0:2,0:3,0:2,0:7),
     &           mOp(4), iuvwx(4)
      Logical JfHss(0:3,0:2,0:3,0:2),JfGrd(0:2,0:3)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      IX(i1,i2)=i1*(i1-1)/2+i2
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
*
c     iRout = 150
c     iPrint = nPrint(iRout)
c     Call qEnter('M1Hss')
*
c     If (iPrint.ge.99) Then
c        Write (6,*) ' In M1Hss: nArr=',nArr
c     End If
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      ipArr = nip
      nArray = nArr - nip +1
*
      iIrrep = 0
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(3,A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      If (la.ge.lb) Then
         call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
      iuvwx(1) = nStab(mdc)
      iuvwx(2) = nStab(ndc)
      mOp(1) = nOp(1)
      mOp(2) = nOp(2)
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
*     Modify the density matrix with the prefactor
*
      nDAO = nElem(la) * nElem(lb)
*     Do 300 iDAO = 1, nDAO
*        Do 310 iZeta = 1, nZeta
*           Fact = Two*rkappa(iZeta)*Pi*ZInv(iZeta)
*           DAO(iZeta,iDAO) = Fact * DAO(iZeta,iDAO)
*310     Continue
*300  Continue
c     If (iPrint.ge.99) Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
*


*     Here we go
*
*-----Loop over nuclear centers
*
      kdc = 0
      Do 100 kCnttp = 1, nCnttp
         If (dbsc(kCnttp)%Charge.eq.Zero) Go To 111
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
            Do 102 lDCRT = 0, nDCRT-1
               Call ICopy(nIrrep*16*9,[0],0,JndHss,1)
               Call iCopy(nIrrep*4*3,[0],0,JndGrd,1)
               Call LCopy(144,[.False.],0,jfHss,1)
               Call LCopy(4,[.False.],0,Tr,1)
               Call LCopy(12,[.False.],0,jfGrd,1)
               mOp(3) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
               mOp(4) = mOp(3)
               TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
               TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
               TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coora(1,3),1)
               call dcopy_(3,TC,1,Coora(1,4),1)
               call dcopy_(3,TC,1,Coori(1,3),1)
               call dcopy_(3,TC,1,Coori(1,4),1)
               If (EQ(A,TC).and.EQ(A,RB)) Goto 102
*
*              COPY CNTLR MATRIXES
*
                Do  iAtom = 0, 1
                  Do iCar  = 0, 2
                    JfGrd(iCar,iAtom) = Ifgrd(iCar,iAtom)
                    Do iIrrep=0,nIrrep-1
                      JndGrd(iCar,iAtom,iIrrep)=
     &                   IndGrd(iCar,iAtom,iIrrep)
                    End Do
                    Do  jAtom = 0, 1
                      Do  jCar = 0, 2
                        JfHss(iAtom,iCar,jAtom,jCar) =
     &                    IfHss(iAtom,iCar,jAtom,jCar)
                        Do iIrrep=0,nIrrep-1
                          JndHss(iAtom,iCar,jAtom,jCar,iIrrep) =
     &                     IndHss(iAtom,iCar,jAtom,jCar,iIrrep)
                        End Do
                      End Do
                    End Do
                  End Do
                End Do

*
                Fact = -dbsc(kCnttp)%Charge*DBLE(nStabM) /
     &             DBLE(LmbdT)
*               Call DYaX(nZeta*nDAO,Fact,DAO,1,Array(ipDAO),1)
                iuvwx(3) = nStab(kdc+kCnt)
                iuvwx(4) = nStab(kdc+kCnt)
*
*-----------Derivatives with respect to the operator is computed via the
*           translational invariance.
*
                nnIrrep=nIrrep
                If (sIrrep) nnIrrep=1
                Do 230 iIrrep=0,nnIrrep-1
                 nDisp = IndDsp(kdc+kCnt,iIrrep)
                 Do 220 iCar = 0, 2
                    iComp = 2**iCar
                    If (TF(kdc+kCnt,iIrrep,iComp)) Then
                       nDisp = nDisp + 1
*
*--------------------Reset flags for the basis set centers so that we
*                    will explicitly compute the derivatives with
*                    respect to those centers. Activate flag for the
*                    third center so that its derivative will be comp-
*                    uted by the translational invariance.
*
                      JndGrd(iCar,0,iIrrep) = Abs(JndGrd(iCar,0,iIrrep))
                      JndGrd(iCar,1,iIrrep) = Abs(JndGrd(iCar,1,iIrrep))
                      JndGrd(iCar,2,iIrrep) = -nDisp
                      JfGrd(iCar,0) = .True.
                      JfGrd(iCar,1) = .True.
                      JfGrd(iCar,2) = .False.
                    Else
                      JndGrd(iCar,2,iIrrep) = 0
                    End If
 220             Continue
 230         Continue
*
*            The third center is claculated by translation invarians
*            This requires the 2nd derivatives on the other centers.
*

             Do iCar=0,2
              Do jAtom=0,2
                if (jAtom.eq.2) Then
                  iStop=iCar
                Else
                  iStop=2
                End If
                Do jCar=0,iStop
                 Do iIrrep=0,nIrrep-1
                  If ((JndGrd(iCar,2,iIrrep).ne.0).and.
     &                (JndGrd(jCar,jAtom,iIrrep).ne.0)) Then
                   JndHss(2,iCar,jAtom,jCar,iIrrep)=
     &                 -IX(Max(Abs(JndGrd(iCar,2,iIrrep)),
     &                 Abs(JndGrd(jCar,jAtom,iIrrep))),
     &                 Min(Abs(JndGrd(iCar,2,iIrrep)),
     &                 Abs(JndGrd(jCar,jAtom,iIrrep))))

                 Tr(2)=.true.
                 If (jAtom.eq.2) Then
                  Maxi=Max(iCar,jCar)
                  Mini=Min(iCar,jCar)
                  jfHss(0,Maxi,0,Mini)=.true.
                  jfHss(1,Maxi,1,Mini)=.true.
                  jfHss(1,iCar,0,jCar)=.true.
                  jfHss(1,jCar,0,iCar)=.true.
                 Else
                  Maxi=Max(iCar,jCar)
                  Mini=Min(iCar,jCar)
                  jfHss(jAtom,Maxi,jAtom,Mini)=.true.
                  jfHss(1,iCar,0,jCar)=.true.
                  jfHss(1,jCar,0,iCar)=.true.
                 End If
                 End If
                 End Do
                End Do
               End Do
              End Do
*
               IfG(0)=.true.
               IfG(1)=.true.
               IfG(2)=.false.
               IfG(3)=.false.
               Do iCent=0,1
                 If (EQ(Coori(1,iCent+1),Coori(1,3) ) ) Then
                 IfG(iCent)=.false.
                 Do iCar=0,2
                   jfGrd(iCar,iCent)=.false.
                   Do kCar=0,2
                    Do KCent=0,3
                     jfHss(iCent,iCar,kCent,kCar)=.false.
                     jfHss(kCent,kCar,iCent,iCar)=.false.
                     Do iIrrep=0,nIrrep-1
                      jndHss(iCent,iCar,kCent,kCar,iIrrep)=0
                      jndHss(kCent,kCar,iCent,iCar,iIrrep)=0
                     End Do
                    End Do
                   End Do
                   Do iIrrep=0,nIrrep-1
                      jndGrd(iCar,iCent,iIrrep)=0
                   End Do
                 End Do
                 End If
               End Do
               Call lCopy(12,[.false.],0,jfgrd,1)
*
               call M1Kernel(Final,Hess,nHess,DAO,nDAO,
     &                   iAnga,nRys,nZeta,
     &                   Array(ipA),Array(ipB),Zeta,ZInv,
     &                   rKappa,P,TC,Coori,Coorac,
     &                   Array(ipArr),nArray,
     &                   jfgrd,jndgrd,jfhss,jndhss,
     &                   ifg,tr,mop,iuvwx,
     &                   kCnttp,Fact,loper(1),0)

*
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
*
c     Call qExit('M1Hss')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
      End If
      End
