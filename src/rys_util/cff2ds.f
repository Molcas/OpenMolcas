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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Cff2DS(nabMax,ncdMax,nRys,
     &                  Zeta,ZInv,Eta,EInv,nT,
     &                  Coori,CoorAC,P,Q,
     &                  la,lb,lc,ld,
     &                  U2,PAQP,QCPQ,B10,B00,lac,B01)
************************************************************************
*                                                                      *
* Object: to compute the coefficients in the three terms recurrence    *
*         relation of the 2D-integrals.                                *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March 1990                                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Zeta(nT), ZInv(nT), Eta(nT), EInv(nT),
     &       Coori(3,4), CoorAC(3,2),
     &       P(nT,3), Q(nT,3), U2(nRys,nT),
     &       PAQP(nRys,nT,3), QCPQ(nRys,nT,3),
     &       B10(nRys,nT,3),
     &       B00(nRys,nT,3),
     &       B01(nRys,nT,3)
*     Local arrays
      Logical AeqB, CeqD, EQ
*define _DEBUG_
#ifdef _DEBUG_
      Character*30 Label
#endif
*
      iRout = 14
      iPrint = nPrint(iRout)
*
#ifdef _DEBUG_
      iPrint=99
      If (iPrint.ge.99) Then
         Call RecPrt(' In Cff2Ds: Coori',' ',Coori,3,4)
         Call RecPrt(' In Cff2Ds: U2',' ',U2,nRys,nT)
      End If
#endif
      AeqB = EQ(Coori(1,1),Coori(1,2))
      CeqD = EQ(Coori(1,3),Coori(1,4))
*
      h12 = Half
      If (nabMax.ne.0 .and. ncdMax.ne.0) Then
         Do 10 iT = 1, nT
            Do iRys = 1, nRys
               B00(iRys,iT,1) = h12 * U2(iRys,iT)
               B10(iRys,iT,1) = ( h12 -
     &            h12 * U2(iRys,iT) * Zeta(iT))*ZInv(iT)
               B01(iRys,iT,1) = B10(iRys,iT,1)
            End Do
 10      Continue
      Else If (ncdMax.eq.0 .and. nabMax.ne.0 .and. lac.eq.0) Then
         Call WarningMessage(2,
     &        'Cff2DS: ncdMax.eq.0 .and. nabMax.ne.0 .and. lac.eq.0')
         Write (6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
         Call Abend()
      Else If (nabMax.eq.0 .and. ncdMax.ne.0 .and. lac.eq.0) Then
         Call WarningMessage(2,
     &        'Cff2DS: nabMax.eq.0 .and. ncdMax.ne.0 .and. lac.eq.0')
         Write (6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
         Call Abend()
      Else If (ncdMax.eq.0 .and. nabMax.ne.0) Then
         Call WarningMessage(2,
     &               'Cff2DS: ncdMax.eq.0 .and. nabMax.ne.0')
         Write (6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
         Call Abend()
      Else If (nabMax.eq.0 .and. ncdMax.ne.0) Then
         Call WarningMessage(2,
     &               'Cff2DS: nabMax.eq.0 .and. ncdMax.ne.0')
         Write (6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
         Call Abend()
      Else  If (nabMax.eq.0 .and. ncdMax.eq.0 .and. lac.ne.0) Then
         Call DYaX(nRys*nT,h12,U2(1,1),1,B00(1,1,1),1)
      End If
      If (nabMax.ne.0) Then
         call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,2),1)
         call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,3),1)
      End If
      If (lac.ne.0) Then
         call dcopy_(nRys*nT,B00(1,1,1),1,B00(1,1,2),1)
         call dcopy_(nRys*nT,B00(1,1,1),1,B00(1,1,3),1)
      End If
      If (ncdMax.ne.0) Then
         call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,2),1)
         call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,3),1)
      End If
*
      If (la+lb.ne.0 .and. lc+ld.ne.0) Then
         If (.Not.AeqB .and. .Not.CeqD) Then
         Do 100 iCar = 1, 3
            Do 110 iT = 1, nT
               Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) = P(iT,iCar) - CoorAC(iCar,1)
                  QCPQ(iRys,iT,iCar) = PAQP(iRys,iT,iCar)
               End Do
 110        Continue
 100     Continue
         Else If (AeqB .and. .Not.CeqD) Then
            Call WarningMessage(2,'Cff2DS: AeqB .and. .Not.CeqD')
            Write (6,*) 'AeqB,CeqD=',AeqB,CeqD
            Call Abend()
         Else If (.Not.AeqB .and. CeqD) Then
            Call WarningMessage(2,'Cff2DS: .Not.AeqB .and. CeqD')
            Write (6,*) 'AeqB,CeqD=',AeqB,CeqD
            Call Abend()
         Else
            call dcopy_(3*nRys*nT,[Zero],0,PAQP,1)
            call dcopy_(3*nRys*nT,[Zero],0,QCPQ,1)
         End If
      Else If (la+lb.ne.0) Then
         Call WarningMessage(2,'Cff2DS: la+lb.ne.0')
         Write (6,*) 'la,lb=',la,lb
         Call Abend()
      Else If (lc+ld.ne.0) Then
         Call WarningMessage(2,'Cff2DS: lc+ld.ne.0')
         Write (6,*) 'lc,ld=',lc,ld
         Call Abend()
      End If
#ifdef _DEBUG_
      If (la+lb.gt.0) Then
         Write (Label,'(A)') ' PAQP(x)'
         Call RecPrt(Label,' ',PAQP(1,1,1),nRys,nT)
         Write (Label,'(A)') ' PAQP(y)'
         Call RecPrt(Label,' ',PAQP(1,1,2),nRys,nT)
         Write (Label,'(A)') ' PAQP(z)'
         Call RecPrt(Label,' ',PAQP(1,1,3),nRys,nT)
      End If
      If (lc+ld.gt.0) Then
         Write (Label,'(A)') ' QCPQ(x)'
         Call RecPrt(Label,' ',QCPQ(1,1,1),nRys,nT)
         Write (Label,'(A)') ' QCPQ(y)'
         Call RecPrt(Label,' ',QCPQ(1,1,2),nRys,nT)
         Write (Label,'(A)') ' QCPQ(z)'
         Call RecPrt(Label,' ',QCPQ(1,1,3),nRys,nT)
      End If
      If (nabMax.ne.0) Then
         Write (Label,'(A)') ' B10(x)'
         Call RecPrt(Label,' ',B10(1,1,1),nRys,nT)
         Write (Label,'(A)') ' B10(y)'
         Call RecPrt(Label,' ',B10(1,1,2),nRys,nT)
         Write (Label,'(A)') ' B10(z)'
         Call RecPrt(Label,' ',B10(1,1,3),nRys,nT)
      End If
      If (lac.ne.0) Then
         Write (Label,'(A)') ' B00(x)'
         Call RecPrt(Label,' ',B00(1,1,1),nRys,nT)
         Write (Label,'(A)') ' B00(y)'
         Call RecPrt(Label,' ',B00(1,1,2),nRys,nT)
         Write (Label,'(A)') ' B00(z)'
         Call RecPrt(Label,' ',B00(1,1,3),nRys,nT)
      End If
      If (ncdMax.ne.0) Then
         Write (Label,'(A)') ' B01(x)'
         Call RecPrt(Label,' ',B01(1,1,1),nRys,nT)
         Write (Label,'(A)') ' B01(y)'
         Call RecPrt(Label,' ',B01(1,1,2),nRys,nT)
         Write (Label,'(A)') ' B01(z)'
         Call RecPrt(Label,' ',B01(1,1,3),nRys,nT)
      End If
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(EInv)
         Call Unused_real_array(Eta)
         Call Unused_real_array(Q)
      End If
      End
