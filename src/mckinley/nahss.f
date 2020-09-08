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
      SubRoutine NAHss(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
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
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, Cff2D
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
     &          iDCRT(0:7),Index(3,4)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nArr), Ccoor(3), Hess(nHess),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
       Logical IfHss(0:1,0:2,0:1,0:2),IfGrd(0:2,0:1), TstFnc, TF,
     &         EQ,IfG(0:3),Tr(0:3)
#ifdef _PATHSCALE_
      Save Fact
#endif
*
*     Local arrrays
*
      Real*8 Coori(3,4), CoorAC(3,2), C(3), TC(3)
      Integer iAnga(4), JndGrd(0:2,0:3,0:7),
     &        JndHss(0:3,0:2,0:3,0:2,0:7),
     &        mOp(4), iuvwx(4)
      Logical JfHss(0:3,0:2,0:3,0:2),JfGrd(0:2,0:3)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      itri(i1,i2)=MAX(i1,i2)*(MAX(i1,i2)-1)/2+MIN(i1,i2)
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,dc(mdc)%nStab)
*
c     iRout = 150
c     iPrint = nPrint(iRout)
c     Call qEnter('NAHSS')
*
c     If (iPrint.ge.99) Then
c        Write (6,*) ' In NAHss: nArr=',nArr
c     End If
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      ipDAO = nip
      nip = nip + nAlpha*nBeta*nElem(la)*nElem(lb)
      If (nip-1.gt.nArr) Then
         Write (6,*) 'NAHss: nip-1.gt.nArr'
         Write (6,*) 'nip,nArr=',nip,nArr
         Call QTrace
         Call Abend()
      End If
      ipArr = nip
      nArray = nArr - nip +1
*
      iIrrep = 0
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
      call dcopy_(3,A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      If (la.ge.lb) Then
         call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
      iuvwx(1) = dc(mdc)%nStab
      iuvwx(2) = dc(ndc)%nStab
      mOp(1) = nOp(1)
      mOp(2) = nOp(2)
*
      ipAOff = ipA
      Do iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
         ipAOff = ipAOff + nAlpha
      End Do
*
      ipBOff = ipB
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
         ipBOff = ipBOff + 1
      End Do
*
*     Modify the density matrix with the prefactor
*
      nDAO = nElem(la) * nElem(lb)
      Do iDAO = 1, nDAO
         Do iZeta = 1, nZeta
            Fact = Two*rkappa(iZeta)*Pi/Zeta(iZeta)
            DAO(iZeta,iDAO) = Fact * DAO(iZeta,iDAO)
         End Do
      End Do
c     If (iPrint.ge.99) Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
*
*-----Loop over nuclear centers
*
      kdc = 0
      Do kCnttp = 1, nCnttp
         If (dbsc(kCnttp)%Charge.eq.Zero) Go To 111
         Do kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)

            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
            Fact = -dbsc(kCnttp)%Charge*DBLE(nStabM) / DBLE(LmbdT)
*
            Call DYaX(nZeta*nDAO,Fact,DAO,1,Array(ipDAO),1)
*
            iuvwx(3) = dc(kdc+kCnt)%nStab
            iuvwx(4) = dc(kdc+kCnt)%nStab
*
            Do 102 lDCRT = 0, nDCRT-1
*
               mOp(3) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
               mOp(4) = mOp(3)
               Call OA(iDCRT(lDCRT),C,TC)
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coori(1,3),1)
               call dcopy_(3,TC,1,Coori(1,4),1)
               If (EQ(A,TC).and.EQ(A,RB)) Goto 102
*
*              Initialize JfGrd, JndGrd, JfHss, and JndHss.
*
               Call LCopy(12,[.False.],0,JfGrd,1)
               Call ICopy(nIrrep*4*3,[0],0,JndGrd,1)
               Call LCopy(144,[.False.],0,JfHss,1)
               Call ICopy(nIrrep*16*9,[0],0,JndHss,1)
*
*              Overwrite with information in IfGrd, IndGrd, IfHss,
*              and IndHss.

               Do iAtom = 0, 1
                  Do iCar  = 0, 2
                     JfGrd(iCar,iAtom) = Ifgrd(iCar,iAtom)
                     Do iIrrep=0,nIrrep-1
                        JndGrd(iCar,iAtom,iIrrep)=
     &                     IndGrd(iCar,iAtom,iIrrep)
                     End Do
                     Do jAtom = 0, 1
                        Do jCar = 0, 2
                           JfHss(iAtom,iCar,jAtom,jCar) =
     &                       IfHss(iAtom,iCar,jAtom,jCar)
                           Do iIrrep=0,nIrrep-1
                              JndHss(iAtom,iCar,jAtom,jCar,iIrrep) =
     &                          IndHss(iAtom,iCar,jAtom,jCar,iIrrep)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
*
*--------------Derivatives with respect to the operator is computed via
*              the translational invariance.
*
               nnIrrep=nIrrep
               If (sIrrep) nnIrrep=1
               Do iIrrep=0,nnIrrep-1
                  nDisp = IndDsp(kdc+kCnt,iIrrep)
                  Do iCar = 0, 2
                     iComp = 2**iCar
                     If (TF(kdc+kCnt,iIrrep,iComp)) Then
                        nDisp = nDisp + 1
*
*-----------------------Reset flags for the basis set centers so that we
*                       will explicitly compute the derivatives with
*                       respect to those centers. Activate flag for the
*                       third center so that its derivative will be comp-
*                       uted by the translational invariance.
*
                        JndGrd(iCar,0,iIrrep)=Abs(JndGrd(iCar,0,iIrrep))
                        JndGrd(iCar,1,iIrrep)=Abs(JndGrd(iCar,1,iIrrep))
                        JndGrd(iCar,2,iIrrep)=-nDisp
                        JfGrd(iCar,0) = .True.
                        JfGrd(iCar,1) = .True.
                        JfGrd(iCar,2) = .False.
                     Else
                        JndGrd(iCar,2,iIrrep) = 0
                     End If
                  End Do
               End Do
*
*              The third center is calculated by translational invariance.
*              This requires the 2nd derivatives on the other centers.
*
               Call LCopy(4,[.False.],0,Tr,1)
               Do iCar=0,2
                  Do jAtom=0,2
                     If (jAtom.eq.2) Then
                        iStop=iCar
                     Else
                        iStop=2
                     End If
                     Do jCar=0,iStop
                        Do iIrrep=0,nIrrep-1
                           If ((JndGrd(iCar,2,iIrrep).ne.0) .and.
     &                         (JndGrd(jCar,jAtom,iIrrep).ne.0)) Then
                              JndHss(2,iCar,jAtom,jCar,iIrrep)=
     &                          -itri(Abs(JndGrd(iCar,2,    iIrrep)),
     &                                Abs(JndGrd(jCar,jAtom,iIrrep)))

                              Tr(2)=.True.
                              If (jAtom.eq.2) Then
                                 Maxi=Max(iCar,jCar)
                                 Mini=Min(iCar,jCar)
                                 JfHss(0,Maxi,0,Mini)=.True.
                                 JfHss(1,Maxi,1,Mini)=.True.
                                 JfHss(1,iCar,0,jCar)=.True.
                                 JfHss(1,jCar,0,iCar)=.True.
                              Else
                                 Maxi=Max(iCar,jCar)
                                 Mini=Min(iCar,jCar)
                                 JfHss(jAtom,Maxi,jAtom,Mini)=.True.
                                 JfHss(1,iCar,0,jCar)=.True.
                                 JfHss(1,jCar,0,iCar)=.True.
                              End If
                           End If
                        End Do
                     End Do
                  End Do
               End Do
*
               IfG(0)=.True.
               IfG(1)=.True.
               IfG(2)=.False.
               IfG(3)=.False.
               Do iCent=0,1
                  If (EQ(Coori(1,iCent+1),Coori(1,3) ) ) Then
                     IfG(iCent)=.False.
                     Do iCar=0,2
                        jfGrd(iCar,iCent)=.False.
                        Do kCar=0,2
                           Do KCent=0,3
                              jfHss(iCent,iCar,kCent,kCar)=.False.
                              jfHss(kCent,kCar,iCent,iCar)=.False.
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
               Call lCopy(12,[.False.],0,jfgrd,1)
*
               nFinal = 0
               Call Rysg2(iAnga,nRys,nZeta,
     &                    Array(ipA),Array(ipB),[One],[One],
     &                    Zeta,ZInv,nZeta,[One],[One],1,
     &                    P,nZeta,TC,1,Coori,Coori,CoorAC,
     &                    Array(ipArr),nArray,
     &                    TNAI1,Fake,Cff2D,
     &                    Array(ipDAO),nDAO,Hess,nHess,
     &                    JfGrd,JndGrd,
     &                    JfHss,JndHss,mOp,iuvwx,ifg,
     &                    nFinal,index,.false.,.true.,tr)
*
 102        Continue
         End Do
 111     kdc = kdc + dbsc(kCnttp)%nCntr
      End Do
*
c     Call qExit('NAHSS')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Final)
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
      End If
      End
