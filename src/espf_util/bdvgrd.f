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
      SubRoutine BdVGrd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                  IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                  iStabM,nStabM)
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "espf.fh"
*
      External TNAI1, Fake, XCff2D
#include "print.fh"
#include "disp.fh"
      Integer IndGrd(3,2), kOp(2), lOper(nComp), iStabM(0:nStabM-1),
     &          iDCRT(0:7)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(4,*),
     &       Array(nZeta*nArr), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
      Logical IfGrad(3,2),ESPFexist
      Character*180 Key,Get_Ln
      External Get_Ln
*
*-----Local arrays
*
      Real*8 C(3), TC(3), Coori(3,4), CoorAC(3,2), ZFd(3),TZFd(3)
      Logical NoLoop, JfGrad(3,4)
      Integer iAnga(4), iStb(0:7),
     &          jCoSet(8,8), JndGrd(3,4), lOp(4), iuvwx(4)
      Character ChOper(0:7)*3
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iPrint = 5
      Call qEnter('BdVGrd')
*
*---- Modify the density matrix with the prefactor
*
      nDAO = nElem(la) * nElem(lb)
      Do iDAO = 1, nDAO
         Do iZeta = 1, nZeta
            Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
            DAO(iZeta,iDAO) = Fact * DAO(iZeta,iDAO)
         End Do
      End Do
      If (iPrint.ge.99) Then
         Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
         Call RecPrt('Alpha',' ',Alpha,1,nAlpha)
         Call RecPrt('Beta ',' ',Beta ,1,nBeta )
         Call RecPrt('Zeta ',' ',Zeta ,1,nZeta )
      End If
*
*---- Loop over charges and dipole moments in the external field
*
      iOrdOp = 0
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      ipDAO = nip
      nip = nip + nAlpha*nBeta*nElem(la)*nElem(lb)*nElem(iOrdOp)
      If (nip-1.gt.nZeta*nArr) Then
         Write (6,*) 'nip-1.gt.nZeta*nArr'
         Call ErrTra
         Call Abend
      End If
      nArray = nZeta*nArr - nip + 1
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = iOrdOp
      iAnga(4) = 0
      call dcopy_(3, A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
*
*     Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
       call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
       call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
      iuvwx(1) = dc(mdc)%nStab
      iuvwx(2) = dc(ndc)%nStab
      lOp(1) = kOp(1)
      lOp(2) = kOp(2)
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
      llOper = lOper(1)
*
*     Loop over centers of the grid
*     But how to retrieve the grid ???
*     I just read it in the ESPF file !
      Call F_Inquire('ESPF.DATA',ESPFExist)
      If (.not.ESPFExist) Then
         Write(6,*)' Error ! No ESPF.DATA file found'
         Call Quit_OnUserError()
      EndIf
      IPotFl = IsFreeUnit(1)
      Call Molcas_Open(IPotFl,'ESPF.DATA')
10    Key = Get_Ln(IPotFl)
      If (Key(1:10).eq.'GRID      ') Then
         Call Get_I1(2,nGrdPt)
         Goto 11
      EndIf
      Goto 10
11    Close (IPotFl)
*
      iDum=0
      Do iPnt = 1, nGrdPt
         ZFd(1)=CCoor(4,iPnt)
         NoLoop = ZFd(1).eq.Zero
         If (NoLoop) Go To 111
*------- Pick up the center coordinates
         C(1)=CCoor(1,iPnt)
         C(2)=CCoor(2,iPnt)
         C(3)=CCoor(3,iPnt)

         If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
*
*------- Generate stabilizor of C
*
         If (nIrrep.eq.8) Then
             nOper=3
         Else If (nIrrep.eq.4) Then
             nOper=2
         Else If (nIrrep.eq.2) Then
             nOper=1
         Else
             nOper=0
         End If
         iChxyz=iChAtm(C,iOper,nOper,iChBas(2))
         Call Stblz(iChxyz,iOper,nIrrep,nStb,iStb,iDum,jCoSet)
*
*--------Find the DCR for M and S
*
         Call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
         Fact = -DBLE(nStabM) / DBLE(LmbdT)
*
         If (iPrint.ge.99) Then
            Write (6,*) ' ZFd=',(ZFd(i),i=1,nElem(iOrdOp))
            Write (6,*) ' Fact=',Fact
            Call RecPrt('DAO*Fact*ZFd()',' ',Array(ipDAO),nZeta*nDAO,
     &                   nElem(iOrdOp))
            Write (6,*) ' m      =',nStabM
            Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),
     &            ii = 0, nStabM-1)
            Write (6,*) ' s      =',nStb
            Write (6,'(9A)') '(S)=',(ChOper(iStb(ii)),
     &            ii = 0, nStb-1)
            Write (6,*) ' LambdaT=',LmbdT
            Write (6,*) ' t      =',nDCRT
            Write (6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),
     &            ii = 0, nDCRT-1)
         End If
         iuvwx(3) = nStb
         iuvwx(4) = nStb
         Call ICopy(6,IndGrd,1,JndGrd,1)
         Do i = 1, 3
            Do j = 1, 2
               JfGrad(i,j) = IfGrad(i,j)
            End Do
         End Do
*
*------- No derivatives with respect to the third or fourth center.
*        The positions of the points in the external field are frozen.
*
         Call ICopy(3,[0],0,JndGrd(1,3),1)
         JfGrad(1,3) = .False.
         JfGrad(2,3) = .False.
         JfGrad(3,3) = .False.
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
         If (iPrint.ge.99) Write (6,*) ' mGrad=',mGrad
         If (mGrad.eq.0) Go To 111
*
         Do lDCRT = 0, nDCRT-1
            lOp(3) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
            lOp(4) = lOp(3)
            Call OA(iDCRT(lDCRT),C,TC)
            call dcopy_(3,TC,1,CoorAC(1,2),1)
            call dcopy_(3,TC,1,Coori(1,3),1)
            call dcopy_(3,TC,1,Coori(1,4),1)
*
            If (iOrdOp.eq.0) Then
               Call DYaX(nZeta*nDAO,Fact*ZFd(1),DAO,1,Array(ipDAO),1)
            Else
               Call OA(iDCRT(lDCRT),ZFd,TZFd)
               jpDAO = ipDAO
               ZFdx=TZFd(1)
               Call DYaX(nZeta*nDAO,Fact*ZFdx,DAO,1,Array(jpDAO),1)
               jpDAO = jpDAO + nZeta*nDAO
               ZFdy=TZFd(2)
               Call DYaX(nZeta*nDAO,Fact*ZFdy,DAO,1,Array(jpDAO),1)
               jpDAO = jpDAO + nZeta*nDAO
               ZFdz=TZFd(3)
               Call DYaX(nZeta*nDAO,Fact*ZFdz,DAO,1,Array(jpDAO),1)
            End If
*
*           Compute integrals with the Rys quadrature.
*
            nT = nZeta
            nDiff=1
            mRys=(la+lb+2+nDiff+iOrdOp)/2
            Call Rysg1(iAnga,mRys,nT,
     &                 Array(ipA),Array(ipB),[One],[One],
     &                 Zeta,ZInv,nZeta,
     &                 [One],[One],1,
     &                 P,nZeta,TC,1,Coori,Coori,CoorAC,
     &                 Array(nip),nArray,
     &                 TNAI1,Fake,XCff2D,
     &                 Array(ipDAO),nDAO*nElem(iOrdOp),
     &                 Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)
*
*           Call RecPrt(' In XFdGrd:Grad',' ',Grad,nGrad,1)
         End Do  ! End loop over DCRs
*
111      Continue
      End Do     ! End loop over centers in the grid
*
      Call qExit('BdVGrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Final)
        Call Unused_integer(nRys)
        Call Unused_integer(nOrdOp)
      End If
      End
