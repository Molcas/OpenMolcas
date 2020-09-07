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
      SubRoutine PPInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: to compute pseudo potential integrals.                       *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "oneswi.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), C(3),
     &       Array(nZeta*nArr), Ccoor(3), TC(3)
      Integer lOper(nComp), iStabM(0:nStabM-1),
     &          iDCRT(0:7), iChO(nComp)
*
      parameter(lproju=9)
      parameter (imax=100,kcrs=1)
      Real*8 ccr(imax),zcr(imax)
      Integer nkcrl(lproju+1,kcrs),nkcru(lproju+1,kcrs),lcr(kcrs),
     &        ncr(imax)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(i) = (i+1)*(i+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*                                                                      *
************************************************************************
*                                                                      *
      nArray=0
      ipScr = 1
      intmax=Max(nElem(la),nElem(lb))
      intmax=intmax**2
      nArray=nArray+intmax
      ipA = ipScr + 2*intmax
      nArray=nArray+nZeta*nElem(la)*nElem(lb)
      If (nArray.gt.nZeta*nArr) Then
         Write (6,*) 'nArray.gt.nZeta*nArr'
         Call QTrace()
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      kdc=-dbsc(1)%nCntr
      Do iCnttp = 1, nCnttp
         kdc = kdc + dbsc(iCnttp)%nCntr

         If (dbsc(iCnttp)%nPP.eq.0) Cycle
cAOM< Get the "true" (non SO) shells
         nPP_S=0
         do kSh = dbsc(iCnttp)%iPP,
     &            dbsc(iCnttp)%iPP + dbsc(iCnttp)%nPP-1
*           Skip if a cardholder shell
            If (Shells(kSh)%nExp.le.0) Cycle
            ncrr=Int(Shells(kSh)%Exp(1))
            if(ncrr.le.500) nPP_S=nPP_S+1
         enddo
         If (nPP_S.eq.0) Cycle
cAOM>
*
         npot = 0
         kShStr=dbsc(iCnttp)%iPP
         kShEnd = kShStr + nPP_S-1
         If (nPP_S-1.gt.lproju) Then
            Write (6,*) 'dbsc(iCnttp)%nPP-1.gt.lproju'
            Write (6,*) 'dbsc(iCnttp)%nPP=',nPP_S
            Write (6,*) 'lproju            =',lproju
            Call QTrace()
            Call Abend()
         End If
         lcr(kcrs)=nPP_S-1
         iSh=0
         iOff = 1
         Do kSh = kShStr, kShEnd
            iSh = iSh + 1
            nkcrl(iSh,kcrs)=iOff
            nkcru(iSh,kcrs)=iOff+Shells(kSh)%nExp/3-1
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
         Do iCntr = 1, dbsc(iCnttp)%nCntr
            C(1:3) = dbsc(iCnttp)%Coor(1:3,iCntr)
*
*
*-----------Find the DCR for M and S
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+iCntr) ,nStab(kdc+iCntr),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            Do lDCRT = 0, nDCRT-1
               Call OA(iDCRT(lDCRT),C,TC)
*                                                                      *
************************************************************************
*                                                                      *
               iZeta = 0
               Do iBeta = 1, nBeta
                  Do iAlpha = 1, nAlpha
                     iZeta = iZeta + 1
                     Call FZero(Array(ipScr),intmax)
                     Call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,
     &                           Beta(iBeta),RB(1),RB(2),RB(3),lb+1,
     &                           Array(ipScr),intmax,Max(la+1,lb+1),
     &                           ccr,zcr,nkcrl,nkcru,lcr,ncr,
     &                           TC(1),TC(2),TC(3),npot)
*
                     Do iB = 1, nElem(lb)
                        Do iA = 1, nElem(la)
                           iAB = (iB-1)*nElem(la)+iA
                           iOff2 = (iB-1)*nElem(la)*nZeta
     &                           + (iA-1)*nZeta + iZeta + ipA - 1
                           Array(iOff2) = Array(iAB+ipScr-1)
                        End Do  ! iA
                     End Do     ! iB
*
                  End Do        ! iAlpha
               End Do           ! iBeta
*                                                                      *
************************************************************************
*                                                                      *
*              Symmetry Adapt
*
               nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
               Call SymAdO(Array(ipA),nZeta,la,lb,nComp,Final,nIC,
     &                     nOp,lOper,iChO,Fact)
            End Do        ! lDCRT
*                                                                      *
************************************************************************
*                                                                      *
         End Do           ! iCntr
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
      End If
      End
