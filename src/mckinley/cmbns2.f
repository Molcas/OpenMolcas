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
* Copyright (C) 1994, Anders Bernhardsson                              *
*               1994, Roland Lindh                                     *
************************************************************************
      SubRoutine CmbnS2(Rnxyz,nZeta,la,lb,Zeta,rKappa,Final,Alpha,Beta,
     &                  Hess,nHess,DAO,IfHss,IndHss,indgrd,iu,iv,nOp)
************************************************************************
*                                                                      *
* Object: compute the 2nd derivative  of the overlap matrix.           *
*                                                                      *
* Called from: OvrHss                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              DDot_   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChTbl
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"

      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), rKappa(nZeta), Beta(nZeta),
     &       Rnxyz(nZeta,3,0:la+2,0:lb+2), Alpha(nZeta),
     &       Hess(nHess),
     &       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
      Logical IfHss(0:1,0:2,0:1,0:2)
      Integer IndHss(0:1,0:2,0:1,0:2,0:nIrrep-1),istb(0:1),
     &          nOp(2), ia(3),
     &          ib(3),indgrd(0:2,0:1,0:nirrep-1)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
*     Index in the triang. local hessian
*
      I(i1,i2)=i1*(i1-1)/2+i2
*
*EAW 970912     ixyz=idLoc(DAO(1,1,1))
c     iRout = 134
c     iPrint = nPrint(iRout)
c     Call qEnter('CmbnS2')
      iStb(0)=iu
      iStb(1)=iv
      iQ = 1
*     Call GetMem(' Enter CmbnS2','LIST','REAL',iDum,iDum)
*
      exp32 = -Three/Two
      Do 25 iZeta = 1, nZeta
         rKappa(iZeta) = rKappa(iZeta) * Zeta(iZeta)**exp32
 25   Continue
c     If (iPrint.ge.99) Then
c        Call RecPrt(' In CmbnS2: Zeta  ',' ',Zeta  ,1,nZeta)
c        Call RecPrt(' In CmbnS2: rKappa',' ',rKappa,1,nZeta)
c        Call RecPrt(' In CmbnS2: Alpha ',' ',Alpha ,1,nZeta)
c        Call RecPrt(' In CmbnS2: Beta  ',' ',Beta  ,1,nZeta)
c     End If
      Do 10 iax = 0, la
         ia(1)=iax
         iyaMax=la-ia(1)
      Do 11 ibx = 0, lb
         ib(1)=ibx
         iybMax=lb-ib(1)
         Do 20 iay = 0, iyaMax
            ia(2)=iay
            ia(3) = la-ia(2)-ia(1)
            ipa= Ind(la,ia(1),ia(3))
         Do 21 iby = 0, iybMax
            ib(2)=iby
          ib(3) = lb-ib(2)-ib(1)
          ipb= Ind(lb,ib(1),ib(3))
*
*
*           Combine overlap integrals
*
*           Integrals like dI/dx1dx1
*
          Do 5 iCoor=1,3
            jCoor=Mod(iCoor,3)+1
            kCoor=Mod(jCoor,3)+1
            If (IfHss(0,iCoor-1,0,iCoor-1)) Then
                  Do 30 iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,I(iCoor,iCoor))=rKappa(iZeta)*
     &                   ((Two*Alpha(iZeta))**2 *
     &                   Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor))*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))-
     &                   Two *  Alpha(iZeta) *
     &                   Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                     If (ia(iCoor).gt.0) Then
                        Final(iZeta,ipa,ipb,I(iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,I(iCoor,iCoor))
     &                      - rKappa(iZeta)*
     &                      (Four *  Alpha(iZeta)* Dble(ia(iCoor)) *
     &                      Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))*
     &                      Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                      Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                     End If
                     If (ia(iCoor).gt.1) Then
                        Final(iZeta,ipa,ipb,I(iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,I(iCoor,iCoor))
     &                      + rKappa(iZeta)*
     &                      (Dble(ia(iCoor)*(ia(iCoor)-1))*
     &                      Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))*
     &                      Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor))*
     &                      Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)))
                     End If
 30               Continue
            End If
 5       Continue
*

*
*           Integrals like dI/dxdz
*
            Do 56 iCoor=2,3
            Do 52 jCoor=1,iCoor-1
                 If (IfHss(0,iCoor-1,0,jCoor-1)) Then
                  Do 51 kCoor=1,3
                    Do 50 iZeta = 1, nZeta
                     If (kCoor.eq.1) Then
                        Final(iZeta,ipa,ipb,
     &                       I(iCoor,jCoor))= rKappa(iZeta)
                     End If
                     If ((kCoor.eq.iCoor).or.(kCoor.eq.jCoor)) Then
                        rIc=Two*Alpha(iZeta)*
     &                      Rnxyz(iZeta,kCoor,ia(kCoor)+1,ib(kCoor))
*
                        If (ia(kCoor).gt.0)
     &                    rIc=rIc-Dble(ia(kCoor))*
     &                       Rnxyz(iZeta,kCoor,ia(kCoor)-1,ib(kCoor))
*
                        Final(iZeta,ipa,ipb,
     &                        I(iCoor,jCoor))=
     &                  Final(iZeta,ipa,ipb,
     &                        I(iCoor,jCoor))*
     &                                       rIc
                     Else
                        Final(iZeta,ipa,ipb,
     &                        I(iCoor,jCoor))=
     &                  Final(iZeta,ipa,ipb,
     &                        I(iCoor,jCoor))*
     &                                Rnxyz(iZeta,kCoor,ia(kCoor),
     &                                ib(kCoor))
                     End If
 50               Continue
 51             Continue
                End If
 52         Continue
 56         Continue
 21      Continue
 20      Continue
 11   Continue
 10   Continue

*
*     Trace the Hessian integrals
*
      nDAO = nZeta * (la+1)*(la+2)/2 * (lb+1)*(lb+2)/2
c     If (iPrint.ge.99) Then
c        Call RecPrt(' S(1)',' ',Final,nDAO,21)
c        Call RecPrt('   D ','(6f12.6)',DAO(1,1,1),nDAO,1)
c     End If
      Do 90 iIrrep=0,nIrrep-1
      Do 100 iCnt=0,1
        Do 105  iCar=1,3
          Do 110 jCnt=0,1
            if (iCnt.eq.jCnt) Then
              iStop=iCar
            Else
              iStop=3
            End If
            Do 115 jCar=1,iStop
            If (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep).ne.0) Then
*
*              Accumulate contribution to the Hessian
*
*
*            Get the characteristics of the diff operator
*
               iCh=iEOr(2**(iCar-1)*iCnt,2**(jCar-1)*jCnt)
*
*             Get the character of the operator in the present irrep
*
               ps=DBLE(iChTbl(iIrrep,nOp(2))**
     &             (iCnt+jCnt))
*
*              Get the transf. character of the diff. operator
*
               ps = ps*DBLE(iPrmt(nOp(2),iCh))
*
*        If the over triangular diff. are needed  multiply by two instead
*        Because of that x2x1 y2y1 z2z1 just appear ones in the (1,2)
*        "subhessian".
*
               if (((iCnt.ne.jCnt).and.(iCar.eq.jCar)).and.
     &               (Abs(indgrd(iCar-1,iCnt,iIrrep)).eq.
     &                Abs(indgrd(jCar-1,jCnt,iIrrep))) ) Then
                         ps=ps*Two
               End If
               iHess = Abs(IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep))
               Fact = DBLE(iStb(iCnt)*iStb(jCnt))/DBLE(nIrrep**2)
               Fact = Fact * ps
               if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep).gt.0) Then
                rtemp=DDot_(nDAO,DAO,1,Final(1,1,1,I(Max(iCar,jCar),
     &                 Min(iCar,jCar))),1)
                Hess(iHess) = Hess(iHess) + Fact*rtemp
               Else
                Fact=Fact*DBLE((-1)**(icnt+jcnt))
                rtemp=DDot_(nDAO,DAO,1,Final(1,1,1,I(Max(iCar,jCar),
     &                 Min(iCar,jCar))),1)
                Hess(iHess) = Hess(iHess) + Fact*rtemp
               End If
            End If
 115       Continue
 110     Continue
 105    Continue
 100  Continue
 90   Continue
*
*     Call GetMem(' Exit CmbnS2','LIST','REAL',iDum,iDum)
c     Call qExit('CmbnS2')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Beta)
      End
