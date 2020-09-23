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
      SubRoutine NAGrd_mck(
#define _CALLING_
#include "grd_mck_interface.fh"
     &                    )
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
      use Center_Info
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, Cff2D
#include "Molcas.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "disp.fh"
#include "disp2.fh"

#include "grd_mck_interface.fh"

*     Local variables
      Integer iDCRT(0:7),Index(3,4)
      Real*8 C(3), TC(3)
      Logical DiffCnt, EQ, Tr(4)
      Real*8 Coora(3,4), Coori(3,4), CoorAC(3,2)
      Integer iAnga(4), JndGrd(3,4,0:7), mOp(4), iuvwx(4),
     &        JndHss(4,3,4,3,0:7), kndgrd(3,4,0:7)
      Logical JfGrd(3,4),kfgrd(3,4),jfg(4), JfHss(4,3,4,3)
      Integer, Parameter:: nPAO=1
      Real*8 :: PAO(nPAO)   ! Dummy array
      Integer, Parameter:: nHess=1
      Real*8 :: Hess(nHess) ! Dummy array
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
c     iRout = 150
c     iPrint = nPrint(iRout)
c     Call qEnter('NAGrd')
*
c     If (iPrint.ge.99) Then
c        Write (*,*) ' In NAGrd: nArr=',nArr
c     End If
      nGrad=lDisp(0)
      Call GetMem('Grad','Allo','REAL',ipGrad,nGrad)
*
      nRys=nHer
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
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
*  Dummies
      Call ICopy(144*nIrrep,[0],0,JndHss,1)
      Call LCopy(144,[.false.],0,jfHss,1)
*
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
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
*-----Loop over nuclear centers
*
      nb=nZeta*nElem(la)*nElem(lb)
      kdc = 0
      Do 100 kCnttp = 1, nCnttp
         If (dbsc(kCnttp)%Charge.eq.Zero) Go To 111
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
            DiffCnt=(IfGrad(iDCar,1).or.IfGrad(iDCar,2))
            If ((.not.DiffCnt).and.((kdc+kCnt).ne.iDCnt)) Goto 101
*
            Call DCR(LmbdT,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
*           Fact = -dbsc(kCnttp)%Charge*DBLE(nStabM*nIrrep) /
*    &             DBLE(LmbdT*dc(kdc+kCnt)%nStab)
            Fact = -dbsc(kCnttp)%Charge*DBLE(nStabM) /
     &             DBLE(LmbdT)
c           If (iPrint.ge.99) Then
c              Write (*,*) ' Charge=',dbsc(kCnttp)%Charge
c              write(*,*)   'NZeta=',nzeta
c              Write(*,*)    'NrOp=',nrop
c              Write (*,*) ' Fact=',Fact
c           End If
            iuvwx(3) = dc(kdc+kCnt)%nStab
            iuvwx(4) = dc(kdc+kCnt)%nStab
            Call LCopy(12,[.false.],0,JFgrd,1)
            Call ICopy(12*nIrrep,[0],0,jndGrd,1)
            Do iCnt = 1, 2
               JfGrd(iDCar,iCnt) = IfGrad(iDCar,iCnt)
            End Do
            Do ICnt=1,2
               If (IfGrad(idcar,iCnt)) Then
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
               mOp(3) = NrOpr(iDCRT(lDCRT))
               mOp(4) = mOp(3)
               Call OA(iDCRT(lDCRT),C,TC)
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coora(1,3),1)
               call dcopy_(3,TC,1,Coora(1,4),1)
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
               Call Rysg2(iAnga,nRys,nZeta,
     &                  Array(ipA),Array(ipB),[One],[One],
     &                  Zeta,ZInv,nZeta,[One],[One],1,
     &                  P,nZeta,TC,1,Coori,Coora,CoorAC,
     &                  Array(nip),nArray,
     &                  TNAI1,Fake,Cff2D,
     &                  PAO,nPAO,Hess,nHess,kfGrd,kndGrd,
     &                  JfHss,JndHss,mOp,iuvwx,Jfg,
     &                  nGr,Index,.true.,.false.,tr)
*
               Do iElem = 1, nElem(la)*nElem(lb)*ngr
                  Do iZeta = 1, nZeta
                     tfac = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
                     indi=(iElem-1)*nZeta+iZeta
                     Array(nip+indi-1) = tfac * Array(nip+indi-1)
                  End Do
               End Do
*
#ifdef _DEBUG_
              Call RecPrt('In NaGrd PI',' ',Array(nip),nb,3)
              Call RecPrt('In NaGrd PI',' ',Final,nb,nrOp)
#endif
               Call SmAdNa(Array(nip),nb,Final,
     &            mop,loper,KndGrd,iuvwx,kfGrd,Index,
     &            idcar,Fact,JFG,tr)
c              IF (iPrint.gt.23)
c    &            Call RecPrt('In NaGrd FI',' ',Final,nb,nrOp)
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
      Call GetMem('Grad','Free','REAL',ipGrad,nGrad)
*
c     Call qExit('NAGrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
         Call Unused_logical_array(Trans)
      End If
      End
