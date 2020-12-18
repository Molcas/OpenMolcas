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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Cff2D(nabMax,ncdMax,nRys,
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
*             March '90                                                *
*                                                                      *
* Modified loop structure for RISC 1991 R. Lindh, Dept. of Theoretical *
* Chemistry, University of Lund, Sweden.                               *
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
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
*     Local arrays
      Character*30 Label
#endif
      Logical AeqB, CeqD, EQ
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' In vCff2D: Coori',' ',Coori,3,4)
      Call RecPrt(' In vCff2D: U2',' ',U2,nRys,nT)
#endif
      AeqB = EQ(Coori(1,1),Coori(1,2))
      CeqD = EQ(Coori(1,3),Coori(1,4))
*
      h12 = Half
      If (nabMax.ne.0 .and. ncdMax.ne.0) Then
         Do iT = 1, nT
            Do iRys = 1, nRys
               B00(iRys,iT,1) = (h12 * U2(iRys,iT))
               B10(iRys,iT,1) = ( h12 -
     &            (h12 * U2(iRys,iT)) *
     &         Eta(iT))*ZInv(iT)
               B01(iRys,iT,1) = ( h12 -
     &            (h12 * U2(iRys,iT)) *
     &         Zeta(iT))*EInv(iT)
            End Do
         End Do
      Else If (ncdMax.eq.0 .and. nabMax.ne.0 .and. lac.eq.0) Then
         Do iT = 1, nT
            Do iRys = 1, nRys
               B10(iRys,iT,1) = ( h12 -
     &            h12 * U2(iRys,iT) * Eta(iT))*ZInv(iT)
            End Do
         End Do
      Else If (nabMax.eq.0 .and. ncdMax.ne.0 .and. lac.eq.0) Then
         Do iT = 1, nT
            Do iRys = 1, nRys
               B01(iRys,iT,1) = ( h12 -
     &            h12 * U2(iRys,iT) * Zeta(iT))*EInv(iT)
            End Do
         End Do
      Else If (ncdMax.eq.0 .and. nabMax.ne.0) Then
         Do iT = 1, nT
            Do iRys = 1, nRys
               B00(iRys,iT,1) = (h12 * U2(iRys,iT))
               B10(iRys,iT,1) = ( h12 -
     &            (h12 * U2(iRys,iT)) * Eta(iT))*ZInv(iT)
            End Do
         End Do
      Else If (nabMax.eq.0 .and. ncdMax.ne.0) Then
         Do iT = 1, nT
            Do iRys = 1, nRys
               B00(iRys,iT,1) = (h12 * U2(iRys,iT))
               B01(iRys,iT,1) = ( h12 -
     &            (h12 * U2(iRys,iT)) * Zeta(iT))*EInv(iT)
            End Do
         End Do
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
            Do iT = 1, nT
               Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) =
     &                P(iT,iCar) - CoorAC(iCar,1) + Eta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
                  QCPQ(iRys,iT,iCar) =
     &                Q(iT,iCar) - CoorAC(iCar,2) - Zeta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               End Do
            End Do
 100     Continue
         Else If (AeqB .and. .Not.CeqD) Then
         Do 200 iCar = 1, 3
            Do iT = 1, nT
               Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) = Eta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
                  QCPQ(iRys,iT,iCar) =
     &                Q(iT,iCar) - CoorAC(iCar,2) - Zeta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               End Do
            End Do
 200     Continue
         Else If (.Not.AeqB .and. CeqD) Then
         Do 300 iCar = 1, 3
            Do iT = 1, nT
               Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) =
     &                P(iT,iCar) - CoorAC(iCar,1) + Eta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
                  QCPQ(iRys,iT,iCar) = - Zeta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               End Do
            End Do
 300     Continue
         Else
         Do 400 iCar = 1, 3
            Do iT = 1, nT
               Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) = Eta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
                  QCPQ(iRys,iT,iCar) = - Zeta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               End Do
            End Do
 400     Continue
         End If
      Else If (la+lb.ne.0) Then
         If (.Not.AeqB) Then
         Do 101 iCar = 1, 3
            Do iT = 1, nT
               Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) = P(iT,iCar) - CoorAC(iCar,1) +
     &                Eta(iT) * U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
               End Do
            End Do
 101     Continue
         Else
         Do 201 iCar = 1, 3
            Do iT = 1, nT
               Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) =
     &                Eta(iT) * U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
               End Do
            End Do
 201     Continue
         End If
      Else If (lc+ld.ne.0) Then
         If (.Not.CeqD) Then
         Do 102 iCar = 1, 3
            Do iT = 1, nT
               Do iRys = 1, nRys
                  QCPQ(iRys,iT,iCar) =
     &                Q(iT,iCar) - CoorAC(iCar,2) +
     &                Zeta(iT) * U2(iRys,iT) * (P(iT,iCar)-Q(iT,iCar))
               End Do
            End Do
 102     Continue
         Else
         Do 202 iCar = 1, 3
            Do iT = 1, nT
               Do iRys = 1, nRys
                  QCPQ(iRys,iT,iCar) =
     &                Zeta(iT) * U2(iRys,iT) * (P(iT,iCar)-Q(iT,iCar))
               End Do
            End Do
 202     Continue
         End If
      End If
#ifdef _DEBUGPRINT_
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
      End
