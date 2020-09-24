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
************************************************************************
      SubRoutine NAGrd(
#define _CALLING_
#include "grd_interface.fh"
     &                )
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
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             October '91                                              *
************************************************************************
      Use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
*     For normal nuclear attraction
      External TNAI1, Fake, Cff2D
*     Finite nuclei
      External TERI1, ModU2, vCff2D
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"

#include "grd_interface.fh"

*     Local variables
      Integer iDCRT(0:7)
      Logical TstFnc, TF
#ifdef _PATHSCALE_
      Save Fact
#endif
      Real*8 Coori(3,4), CoorAC(3,2), C(3), TC(3)
      Integer iAnga(4), JndGrd(3,4), lOp(4), iuvwx(4)
      Logical JfGrad(3,4)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,
     &                       dc(mdc)%nStab)
*
      iRout = 150
      iPrint = nPrint(iRout)
*     Call qEnter('NAGrd')
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Write (6,*) ' In NAGrd: nArr=',nArr
         nDAO = nElem(la) * nElem(lb)
         Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
      End If
#endif
*
      nRys=nHer
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      ipDAO = nip
      nip = nip + nAlpha*nBeta*nElem(la)*nElem(lb)
      If (nip-1.gt.nZeta*nArr) Then
         Write (6,*) ' nip-1.gt.nZeta*nArr'
         Call Abend()
      End If
      nArray = nZeta*nArr - nip +1
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
*     Modify the density matrix with the prefactor
*
      nDAO = nElem(la) * nElem(lb)
      If (Nuclear_Model.eq.Point_Charge) Then
         Do iDAO = 1, nDAO
            Do iZeta = 1, nZeta
               Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
               DAO(iZeta,iDAO) = Fact * DAO(iZeta,iDAO)
            End Do
         End Do
      End If
C     If (iPrint.ge.99) Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
*
*-----Loop over nuclear centers
*
      kdc = 0
      Do kCnttp = 1, nCnttp
         If (dbsc(kCnttp)%Charge.eq.Zero) Go To 111
         Do kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
            Fact = -dbsc(kCnttp)%Charge*DBLE(nStabM) / DBLE(LmbdT)
*
*           Modify the density matrix with prefactors in case of finite nuclei
*
            If (Nuclear_Model.eq.Gaussian_Type) Then
               Eta=dbsc(kCnttp)%ExpNuc
               rKappcd=TwoP54/Eta
*              Tag on the normalization factor of the nuclear Gaussian
               Fact=Fact*(Eta/Pi)**(Three/Two)
               jpDAO=ipDAO
               Do iDAO = 1, nDAO
                  Do iZeta = 1, nZeta
*                    On flight modification of Kappa
                     rKappab=TwoP54*rKappa(iZeta)/Zeta(iZeta)
                     Array(jpDAO)=Fact*DAO(iZeta,iDAO)*rKappab*rKappcd
     &                           *Sqrt(One/(Zeta(iZeta)+Eta))
                     jpDAO=jpDAO+1
                  End Do
               End Do
            Else If (Nuclear_Model.eq.Point_Charge) Then
               Call DYaX(nZeta*nDAO,Fact,DAO,1,Array(ipDAO),1)
            Else
               Write (6,*) 'NaGrd: Fermi type nuclear distribution not '
     &                   //'implemented yet!'
               Call Abend()
            End If
            iuvwx(3) = dc(kdc+kCnt)%nStab
            iuvwx(4) = dc(kdc+kCnt)%nStab
            Call ICopy(6,IndGrd,1,JndGrd,1)
            Do i = 1, 3
               Do j = 1, 2
                  JfGrad(i,j) = IfGrad(i,j)
               End Do
            End Do
*
*-----------Derivatives with respect to the operator is computed via the
*           translational invariance.
*
            nDisp = IndDsp(kdc+kCnt,iIrrep)
            Do iCar = 0, 2
               iComp = 2**iCar
               If ( TF(kdc+kCnt,iIrrep,iComp) .and.
     &              .Not.dbsc(kCnttp)%Frag .and.
     &              .Not.dbsc(kCnttp)%pChrg ) Then
                  nDisp = nDisp + 1
                  If (Direct(nDisp)) Then
*--------------------Reset flags for the basis set centers so that we
*                    will explicitly compute the derivatives with
*                    respect to those centers. Activate flag for the
*                    third center so that its derivative will be comp-
*                    uted by the translational invariance.
                     JndGrd(iCar+1,1) = Abs(JndGrd(iCar+1,1))
                     JndGrd(iCar+1,2) = Abs(JndGrd(iCar+1,2))
                     JndGrd(iCar+1,3) = -nDisp
                     JfGrad(iCar+1,1) = .True.
                     JfGrad(iCar+1,2) = .True.
                     JfGrad(iCar+1,3) = .False.
                  Else
                     JndGrd(iCar+1,3) = 0
                     JfGrad(iCar+1,3) = .False.
                  End If
               Else
                  JndGrd(iCar+1,3) = 0
                  JfGrad(iCar+1,3) = .False.
               End If
            End Do
*-----------No derivatives with respect to the fourth center.
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
C           If (iPrint.ge.99) Write (6,*) ' mGrad=',mGrad
            If (mGrad.eq.0) Go To 101
*
            Do lDCRT = 0, nDCRT-1
               lOp(3) = NrOpr(iDCRT(lDCRT))
               lOp(4) = lOp(3)
               Call OA(iDCRT(lDCRT),C,TC)
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coori(1,3),1)
               call dcopy_(3,TC,1,Coori(1,4),1)
*
*
               If (Nuclear_Model.eq.Gaussian_Type) Then
                  Eta=dbsc(kCnttp)%ExpNuc
                  EInv=One/Eta
                  Call Rysg1(iAnga,nRys,nZeta,
     &                       Array(ipA),Array(ipB),[One],[One],
     &                       Zeta,ZInv,nZeta,[Eta],[EInv],1,
     &                       P,nZeta,TC,1,Coori,Coori,CoorAC,
     &                       Array(nip),nArray,
     &                       TERI1,ModU2,vCff2D,
     &                       Array(ipDAO),nDAO,Grad,nGrad,
     &                       JfGrad,JndGrd,lOp,iuvwx)
               Else If (Nuclear_Model.eq.Point_Charge) Then
                  Call Rysg1(iAnga,nRys,nZeta,
     &                       Array(ipA),Array(ipB),[One],[One],
     &                       Zeta,ZInv,nZeta,[One],[One],1,
     &                       P,nZeta,TC,1,Coori,Coori,CoorAC,
     &                       Array(nip),nArray,
     &                       TNAI1,Fake,Cff2D,
     &                       Array(ipDAO),nDAO,Grad,nGrad,
     &                       JfGrad,JndGrd,lOp,iuvwx)
               Else
*...more to come...
               End If
*
C              Call RecPrt('In NaGrd: Grad',' ',Grad,nGrad,1)
            End Do
 101     Continue
         End Do
 111     kdc = kdc + dbsc(kCnttp)%nCntr
      End Do
*
*     Call qExit('NAGrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Final)
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
      End If
      End
