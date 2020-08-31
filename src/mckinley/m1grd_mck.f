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
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine m1Grd_mck(Alpha,nAlpha,Beta, nBeta,
     &                 Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,la,lb,A,RB,nRys,
     &                 Array,nArr,Ccoor,nOrdOp,
     &                 IfGrd,IndGrd,nOp,
     &                 lOper,iu,iv,nrOp,iDCar,iDCnt,iStabM,nStabM,ldum)
************************************************************************
*                                                                      *
* Object: to compute the gradient of the nuclear attraction integrals. *
*          Something is wrong here                                     *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy  (ESSL)                                           *
*              ICopy                                                   *
*              Rysg1                                                   *
*              QExit                                                   *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             October 1991                                             *
*              Anders Bernhardsson 1995                                *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, Cff2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "disp.fh"
#include "disp2.fh"
      Integer IndGrd(0:nIrrep-1), nOp(2),
     &          iDCRT(0:7),inddum(144*8)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nrOp),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), C(3), TC(3),
     &       Array(nArr),cCoor(3)
      Logical IfGrd(3,2), DiffCnt,EQ,ldum(2),Tr(4),ifdum(144)
*
*     Local arrrays
*
      Real*8 Coori(3,4), CoorAC(3,2)
      Integer iAng(4), JndGrd(3,4,0:7),
     &          mOp(4), iuvwx(4),
     &          kndgrd(3,4,0:7),iStabM(0:7)
      Logical JfGrd(3,4),kfgrd(3,4),jfg(4)
      Dimension Dum(1)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
*
c     If (iPrint.ge.99) Then
c        Write (*,*) ' In NAGrd: nArr=',nArr
c     End If
      nGrad=lDisp(0)
      Call GetMem('Grad','Allo','REAL',ipGrad,nGrad)
      call icopy(144*nirrep,[0],0,inddum,1)
      call lcopy(144,[.false.],0,ifdum,1)
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      If (nip-1.gt.nArr)
     &   Write (6,*) ' nip-1.gt.nArr'
      nArray = nArr - nip +1
*
      iIrrep = 0
      iAng(1) = la
      iAng(2) = lb
      iAng(3) = 0
      iAng(4) = 0
*  Dummies
*
      call dcopy_(3,A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      If (la.ge.lb) Then
         call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
      iuvwx(1) = iu
      iuvwx(2) = iv
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
      nDAO = 1
*
*-----Loop over nuclear centers
*
      nb=nZeta*nElem(la)*nElem(lb)
      kdc = 0
      Do 100 kCnttp = 1, nCnttp
         If (.Not.dbsc(kCnttp)%ECP) Go To 111
         If (dbsc(kCnttp)%nM1.eq.0) Go To 111

         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
            DiffCnt=(IfGrd(iDCar,1).or.IfGrd(iDCar,2))
            If ((.not.DiffCnt).and.((kdc+kCnt).ne.iDCnt)) Goto 101
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
*           Fact = -Charge(kCnttp)*DBLE(nStabM*nIrrep) /
*    &             DBLE(LmbdT*nStab(kdc+kCnt))
            Fact = -Charge(kCnttp)*DBLE(nStabM) /
     &             DBLE(LmbdT)
c           If (iPrint.ge.99) Then
c              Write (*,*) ' Charge=',Charge(kCnttp)
c              write(*,*)   'NZeta=',nzeta
c              Write(*,*)    'NrOp=',nrop
c              Write (*,*) ' Fact=',Fact
c           End If
            iuvwx(3) = nStab(kdc+kCnt)
            iuvwx(4) = nStab(kdc+kCnt)
            Call LCopy(12,[.false.],0,JFgrd,1)
            Call ICopy(12*nIrrep,[0],0,jndGrd,1)
            Do iCnt = 1, 2
                  JfGrd(iDCar,iCnt) = IfGrd(iDCar,iCnt)
            End Do
            Do ICnt=1,2
               If (ifgrd(idcar,iCnt)) Then
                 Do iIrrep=0,nIrrep-1
                   jndGrd(iDCar,iCnt,iIrrep)=IndGrd(iIrrep)
                 End Do
               End IF
            End Do
*
            Tr(1)=.false.
            Tr(2)=.false.
            Tr(3)=.false.
            Tr(4)=.false.
            If ((kdc+kCnt).eq.iDCnt) Then
                 Tr(3)=.true.
                 JfGrd(iDCar,1) = .true.
                 JfGrd(iDCar,2) = .true.
                 Do iIrrep=0,nIrrep-1
                 jndGrd(iDCar,3,iIrrep) = - IndGrd(iIrrep)
                 End Do
            End If
*
            Do 102 lDCRT = 0, nDCRT-1
               Call lCopy(12,JfGrd,1,kfGrd,1)
               Call iCopy(12*nIrrep,JndGrd,1,kndgrd,1)
               mOp(3) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
               mOp(4) = mOp(3)
               TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
               TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
               TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coori(1,3),1)
               call dcopy_(3,TC,1,Coori(1,4),1)
               If (Eq(A,RB).and.EQ(A,TC)) goto 102
               If (EQ(A,TC)) Then
                 kfGrd(iDCar,1) = .false.
                 Do iIrrep=0,nIrrep-1
                  kndGrd(iDCar,1,iirrep)=0
                 End Do
               End If
               If (EQ(RB,TC)) Then
                 kfGrd(iDCar,2) = .false.
                 Do iIrrep=0,nIrrep-1
                  kndgrd(iDCar,2,iIrrep)=0
                 End Do
               End If
*
               If (kfGrd(idcar,1)) Then
                JFG(1)=.true.
               Else
                JFG(1)=.false.
               End If
               If (kfGrd(idcar,2)) Then
                JFG(2)=.true.
               Else
                JFG(2)=.false.
               End If
               JFG(3)=.false.
               JFG(4)=.false.




               call M1Kernel(Final,Dum,0,Dum,0,
     &                   iAng,nRys,nZeta,
     &                   Array(ipA),Array(ipB),Zeta,ZInv,
     &                   rKappa,P,TC,Coori,Coorac,
     &                   Array(nip),nArr-nip+1,
     &                   kfgrd,kndgrd,ifdum,inddum,
     &                   jfg,tr,mop,iuvwx,
     &                   kCnttp,fact,loper,idcar)

 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
      Call GetMem('Grad','Free','REAL',ipGrad,nGrad)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
         Call Unused_logical_array(ldum)
      End If
      End
