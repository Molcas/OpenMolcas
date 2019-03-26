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
* Copyright (C) 1995,2001, Roland Lindh                                *
*               2003, Michael Diedenhofen                              *
************************************************************************
      SubRoutine COSGrd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                  IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                  iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of electronic COSMO cont. *
*         integrals.                                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              DCR                                                     *
*              XRysg1                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*             M. Diedenhofen Nov. 2003                                 *
*             changes pcmgrd routines which do not take into account   *
*             the contribution of a non fixed grid                     *
*     Author of orig routine: Roland Lindh,                            *
*             Dept. of Theoretical Chemistry, University               *
*             of Lund, Sweden, May '95                                 *
*                                                                      *
*             Modified to PCM gradients September 2001, Lund, by       *
*             R. Lindh.                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, Cff2D
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "rctfld.fh"

      Integer IndGrd(3,2), kOp(2), lOper(nComp), iStabM(0:nStabM-1),
     &          iDCRT(0:7)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(3,nComp),
     &       Array(nZeta*nArr), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
      Logical IfGrad(3,2)
*
*-----Local arrys
*
      Real*8 C(3), TC(3), Coori(3,4), CoorAC(3,2)
      Logical NoLoop, JfGrad(3,4)
      Integer iAnga(4), iStb(0:7), JndGrd(3,4), lOp(4), iuvwx(4)
      Character ChOper(0:7)*3
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 151
      iPrint = nPrint(iRout)
      Call qEnter('COSgrd')
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
      If (iPrint.ge.99) Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      ipDAO = nip
      nip = nip + nAlpha*nBeta*nElem(la)*nElem(lb)*nElem(nOrdOp)
      If (nip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2,'Error in COSGrd')
         Write (6,*) 'nip-1.gt.nZeta*nArr'
         Call Abend()
      End If
      nArray = nZeta*nArr - nip + 1
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = nOrdOp
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
      iuvwx(1) = nStab(mdc)
      iuvwx(2) = nStab(ndc)
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
*     Loop over the tiles
*
      Do iTs = 1, nTs
         Q=Work((iTs-1)+ip_Q)
         NoLoop = Q.eq.Zero
         If (NoLoop) Go To 111
*------- Pick up the tile coordinates
         C(1) = Work((iTs-1)*4+ip_Tess  )
         C(2) = Work((iTs-1)*4+ip_Tess+1)
         C(3) = Work((iTs-1)*4+ip_Tess+2)
         kat  = nint(Work((iTs-1)*4+ip_Tess+3))
         If (iPrint.ge.99) Call RecPrt('C',' ',C,3,1)
*
*------- Generate stabilizor of C
*
         nStb=1
         iStb(0)=0
*
*--------Find the DCR for M and S
*
         Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &            iStb,nStb,iDCRT,nDCRT)
         Fact = -DBLE(nStabM) / DBLE(LmbdT)
*
         If (iPrint.ge.99) Then
            Write (6,*) ' Q=',Q
            Write (6,*) ' Fact=',Fact
            Call RecPrt('DAO*Fact*Q',' ',Array(ipDAO),nZeta*nDAO,
     &                   nElem(nOrdOp))
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
c        skip contributions if the segment belongs to the atom
         if(mdc.eq.kat) then
c             skip 1 center
              JfGrad(1,1) = .false.
              JfGrad(2,1) = .false.
              JfGrad(3,1) = .false.
         endif
         if(ndc.eq.kat) then
c             skip 2 center
              JfGrad(1,2) = .false.
              JfGrad(2,2) = .false.
              JfGrad(3,2) = .false.
         endif
*
*         set up third center
*
          iIrrep=0
          nDisp = IndDsp(kat,iIrrep)
          Do iCar = 0, 2
             nDisp = nDisp + 1
             JndGrd(iCar+1,1) = Abs(JndGrd(iCar+1,1))
             JndGrd(iCar+1,2) = Abs(JndGrd(iCar+1,2))
             JndGrd(iCar+1,3) = -nDisp
             JfGrad(iCar+1,3) = .False.
          enddo

         Call ICopy(3,[0],0,JndGrd(1,4),1)
         JfGrad(1,4) = .False.
         JfGrad(2,4) = .False.
         JfGrad(3,4) = .False.

         mGrad = 0
         Do iCar = 1, 3
            Do i = 1, 3
               If (JfGrad(iCar,i)) mGrad = mGrad + 1
            End Do
         End Do
         If (iPrint.ge.99) Write (6,*) ' mGrad=',mGrad
         If (mGrad.eq.0) Go To 111
*
         Do lDCRT = 0, nDCRT-1
            lOp(3) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
            lOp(4) = lOp(3)
            TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
            TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
            TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
            call dcopy_(3,TC,1,CoorAC(1,2),1)
            call dcopy_(3,TC,1,Coori(1,3),1)
            call dcopy_(3,TC,1,Coori(1,4),1)
*
            Call DYaX(nZeta*nDAO,Fact*Q,DAO,1,Array(ipDAO),1)
*
*           Compute integrals with the Rys quadrature.
*
            nT = nZeta
            nDiff=1
            mRys=nRys

            Call Rysg1(iAnga,mRys,nT,
     &                 Array(ipA),Array(ipB),[One],[One],
     &                 Zeta,ZInv,nZeta,
     &                 [One],[One],1,
     &                 P,nZeta,TC,1,Coori,Coori,CoorAC,
     &                 Array(nip),nArray,
     &                 TNAI1,Fake,Cff2D,
     &                 Array(ipDAO),nDAO*nElem(nOrdOp),
     &                 Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)
*
*           Call RecPrt(' In COSgrd:Grad',' ',Grad,nGrad,1)
         End Do  ! End loop over DCRs
*
111      Continue
      End Do     ! End loop over centers in the external field
*
*     Call GetMem(' Exit COSgrd','LIST','REAL',iDum,iDum)
      Call qExit('COSgrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Final)
         Call Unused_real_array(Ccoor)
      End If
      End
