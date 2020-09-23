************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine PPGrd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                  IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                  iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: to compute pseudo potential gradient integrals               *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
************************************************************************
      Use Basis_Info
      use Center_Info
      use Symmetry_Info, only: iOper
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "oneswi.fh"
#include "print.fh"
#include "disp.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), C(3),
     &       Array(nZeta*nArr), Ccoor(3), TC(3),
     &       Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
      Character*80 Label
      Integer lOper(nComp), iStabM(0:nStabM-1),
     &        iDCRT(0:7), IndGrd(3,2), iuvwx(4),
     &        kOp(2), lOp(4), JndGrd(3,4)
      Logical IfGrad(3,2), JfGrad(3,4)
*
      parameter(lproju=9)
      parameter (imax=100,kcrs=1)
      Real*8 ccr(imax),zcr(imax)
      Integer nkcrl(lproju+1,kcrs),nkcru(lproju+1,kcrs),lcr(kcrs),
     &        ncr(imax)
*
      Logical EQ, TstFnc, TF
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(i) = (i+1)*(i+2)/2
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,dc(mdc)%nStab)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 122
      iPrint = nPrint(iRout)
*
      nDAO=nElem(la)*nElem(lb)
      iIrrep = 0
      iuvwx(1) = dc(mdc)%nStab
      iuvwx(2) = dc(ndc)%nStab
      lOp(1) = iOper(kOp(1))
      lOp(2) = iOper(kOp(2))
*                                                                      *
************************************************************************
*                                                                      *
*     Memory managment
*
      nArray=0
      ipRef=1
*
*     la+1, lb
*
      nlaplb=Max(nElem(la+1),nElem(lb))**2
      iplaplb = ipRef + 2*nArray
      nArray=nArray+nlaplb
*
*     la-1, lb
*
      If (la.gt.0) Then
         nlamlb=Max(nElem(la-1),nElem(lb))**2
      Else
         nlamlb=0
      End If
      iplamlb = ipRef + 2*nArray
      nArray=nArray+nlamlb
*
*     la, lb+1
*
      nlalbp=Max(nElem(la),nElem(lb+1))**2
      iplalbp = ipRef + 2*nArray
      nArray=nArray+nlalbp
*
*     la, lb-1
*
      If (lb.gt.0) Then
         nlalbm=Max(nElem(la),nElem(lb-1))**2
      Else
         nlalbm=0
      End If
      iplalbm = ipRef + 2*nArray
      nArray=nArray+nlalbm
*
      If (nArray.gt.nZeta*nArr) Then
         Write (6,*) 'nArray.gt.nZeta*nArr'
         Call QTrace()
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      iComp = 1
      kdc =-dbsc(1)%nCntr
      Do iCnttp = 1, nCnttp
         kdc = kdc + dbsc(iCnttp)%nCntr
         If (dbsc(iCnttp)%Charge.eq.0d0) Cycle
         If (dbsc(iCnttp)%nPP.eq.0) Cycle
cAOM< Get the "true" (non SO) shells
         nPP_S=0
         do kSh = dbsc(iCnttp)%iPP,
     &            dbsc(iCnttp)%iPP + dbsc(iCnttp)%nPP-1
           If (Shells(kSh)%nExp.le.0) Cycle
           ncrr=Int(Shells(kSh)%Exp(1))
           if(ncrr.le.500) nPP_S=nPP_S+1
         enddo
         If (nPP_S.eq.0) Cycle
cAOM>
*
         npot = 0
         kShStr=dbsc(iCnttp)%iPP
cAOM         kShEnd = kShStr + dbsc(iCnttp)%nPP-1
         kShEnd = kShStr + nPP_S-1
         If (dbsc(iCnttp)%nPP-1.gt.lproju) Then
            Write (6,*) 'dbsc(iCnttp)%nPP-1.gt.lproju'
CAOM            Write (6,*) 'dbsc(iCnttp)%nPP=',dbsc(iCnttp)%nPP
            Write (6,*) 'dbsc(iCnttp)%nPP=',nPP_S
            Write (6,*) 'lproju            =',lproju
            Call QTrace()
            Call Abend()
         End If
CAOM         lcr(kcrs)=dbc(iCnttp)%nPP-1
         lcr(kcrs)=nPP_S-1
         iSh=0
         iOff = 1
         Do kSh = kShStr, kShEnd
            iSh = iSh + 1
            nkcrl(iSh,kcrs)=iOff
            nkcru(iSh,kcrs)=iOff+Shells(ksh)%nExp/3-1
            iOff = iOff + Shells(kSh)%nExp/3
            If (nPot.gt.imax) Then
               Write (6,*)' Pseudo: nPot.gt.imax'
               Write (6,*)'         nPot=',nPot
               Write (6,*)'         imax=',imax
               Call QTrace()
               Call Abend()
            End If
            iStrt=1
            Do iExp = 1, Shells(kSh)%nExp/3
               npot = npot + 1
               ncr(npot)=Int(Shells(kSh)%Exp(iStrt  ))
               zcr(npot)=    Shells(kSh)%Exp(iStrt+1)
               ccr(npot)=    Shells(kSh)%Exp(iStrt+2)
               iStrt=iStrt+3
            End Do
         End Do
C        Write (*,*) 'ncr',(ncr(i),i=1,npot)
C        Write (*,*) 'zcr',(zcr(i),i=1,npot)
C        Write (*,*) 'ccr',(ccr(i),i=1,npot)
C        Write (*,*) 'nkcrl',(nkcrl(i,1),i=1,iSh)
C        Write (*,*) 'nkcru',(nkcru(i,1),i=1,iSh)
*
         Do kCnt = 1, dbsc(iCnttp)%nCntr
            C(1:3) = dbsc(iCnttp)%Coor(1:3,kCnt)
*
*-----------Find the DCR for M and S
*
            Call DCR(LmbdT,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab ,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            iuvwx(3) = dc(kdc+kCnt)%nStab
            iuvwx(4) = dc(kdc+kCnt)%nStab
            Call ICopy(6,IndGrd,1,JndGrd,1)
            Do i = 1, 3
               Do j = 1, 2
                  JfGrad(i,j) = IfGrad(i,j)
               End Do
            End Do
*
            nDisp = IndDsp(kdc+kCnt,iIrrep)
            Do iCar = 0, 2
               JfGrad(iCar+1,3) = .False.
               iCmp = 2**iCar
               If ( TF(kdc+kCnt,iIrrep,iCmp) .and.
     &              .Not.dbsc(iCnttp)%pChrg ) Then
                  nDisp = nDisp + 1
                  If (Direct(nDisp)) Then
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
            End Do
            Call ICopy(3,[0],0,JndGrd(1,4),1)
            JfGrad(1,4) = .False.
            JfGrad(2,4) = .False.
            JfGrad(3,4) = .False.
            mGrad = 0
            Do iCar = 1, 3
               Do i = 1, 2
                  If (JfGrad(iCar,i)) mGrad = mGrad + 1
               End Do
            End Do
            If (mGrad.eq.0) Cycle
*
            Do lDCRT = 0, nDCRT-1
               lOp(3) = iDCRT(lDCRT)
               lOp(4) = lOp(3)
               Call OA(iDCRT(lDCRT),C,TC)
               If (EQ(A,RB).and.EQ(A,TC)) Cycle
*                                                                      *
************************************************************************
*                                                                      *
               iZeta = 0
               Do iBeta = 1, nBeta
                  Do iAlpha = 1, nAlpha
                     iZeta = iZeta + 1
*
*                    la+1, lb
*
                     Call FZero(Array(iplaplb),nlaplb)
                     Call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+2,
     &                           Beta(iBeta),RB(1),RB(2),RB(3),lb+1,
     &                           Array(iplaplb),nlaplb,Max(la+2,lb+1),
     &                           ccr,zcr,nkcrl,nkcru,lcr,ncr,
     &                           TC(1),TC(2),TC(3),npot)
*
*                    la-1, lb
*
                     If (la.gt.0) Then
                     Call FZero(Array(iplamlb),nlamlb)
                     Call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la  ,
     &                           Beta(iBeta),RB(1),RB(2),RB(3),lb+1,
     &                           Array(iplamlb),nlamlb,Max(la  ,lb+1),
     &                           ccr,zcr,nkcrl,nkcru,lcr,ncr,
     &                           TC(1),TC(2),TC(3),npot)
                     End If
*
*                    la, lb+1
*
                     Call FZero(Array(iplalbp),nlalbp)
                     Call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,
     &                           Beta(iBeta),RB(1),RB(2),RB(3),lb+2,
     &                           Array(iplalbp),nlalbp,Max(la+1,lb+2),
     &                           ccr,zcr,nkcrl,nkcru,lcr,ncr,
     &                           TC(1),TC(2),TC(3),npot)
*
*                    la, lb-1
*
                     If (lb.gt.0) Then
                     Call FZero(Array(iplalbm),nlalbm)
                     Call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,
     &                           Beta(iBeta),RB(1),RB(2),RB(3),lb  ,
     &                           Array(iplalbm),nlalbm,Max(la+1,lb  ),
     &                           ccr,zcr,nkcrl,nkcru,lcr,ncr,
     &                           TC(1),TC(2),TC(3),npot)
                     End If
*
*                    Assemble gradient and store in Final.
*
                     Call Assemble_PPGrd(Final,nZeta,la,lb,
     &                                   iZeta,
     &                                   Alpha(iAlpha),Beta(iBeta),
     &                                   Array(iplaplb),
     &                                   Array(iplamlb),
     &                                   Array(iplalbp),
     &                                   Array(iplalbm),JfGrad)
*
                  End Do        ! iAlpha
               End Do           ! iBeta
*
CAOM<
               If(Abs(Fact-1d0).gt.1d-7)
     &         call dscal_(nAlpha*nBeta*nElem(la)*nElem(lb)*mGrad,
     &                   Fact,Final,1)
CAOM>
               If (iPrint.ge.99) Then
                  Write (6,*) ' Result in PPGrd'
                  Write (6,*) JfGrad
                  Do ia = 1, nElem(la)
                     Do ib = 1, nElem(lb)
                        Do iVec = 1, mGrad
                           Write (Label,'(A,I2,A,I2,A)')
     &                           ' Final(',ia,',',ib,')'
                           Call RecPrt(Label,' ',Final(1,ia,ib,iVec),
     &                                 nAlpha,nBeta)
                        End Do
                     End Do
                  End Do
               End If
*                                                                      *
************************************************************************
*                                                                      *
*              Distribute contributions to the gradient
*
               Call Distg1X(Final,DAO,nZeta,nDAO,mGrad,Grad,nGrad,
     &                     JfGrad,JndGrd,iuvwx,lOp)
*
            End Do        ! lDCRT
*                                                                      *
************************************************************************
*                                                                      *
         End Do           ! kCnt
      End Do              ! iCnttp
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Zeta)
         Call Unused_real_array(ZInv)
         Call Unused_real_array(rKappa)
         Call Unused_real_array(P)
         Call Unused_integer(nHer)
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
      End If
      End
