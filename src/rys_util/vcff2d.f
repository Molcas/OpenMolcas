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
* Copyright (C) 1990,1991,1994, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine vCff2D(iDum1 ,iDum2 ,nRys,
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
*             Modified loop structure for RISC 1991                    *
*             Modified for decreased memory access January '94.        *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Zeta(nT), ZInv(nT), Eta(nT), EInv(nT),
     &       Coori(3,4), CoorAC(3,2),
     &       P(nT,3), Q(nT,3), U2(nRys,nT),
     &       PAQP(nRys,nT,3), QCPQ(nRys,nT,3),
     &       B10(nRys,nT),
     &       B00(nRys,nT),
     &       B01(nRys,nT)
      Real*8 tmp
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
*     Local arrays
      Character*30 Label
#endif
      Logical AeqB, CeqD, EQ, PrintB10, PrintB00, PrintB01
*
      iRout = 14
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' In vCff2D: Coori',' ',Coori,3,4)
      Call RecPrt(' In vCff2D: U2',' ',U2,nRys,nT)
      Call RecPrt(' in vCff2d: Zeta',' ',Zeta,1,nT)
      Call RecPrt(' in vCff2d: Eta ',' ',Eta, 1,nT)
      Call RecPrt(' in vCff2d: ZInv',' ',ZInv,1,nT)
      Call RecPrt(' in vCff2d: EInv',' ',EInv,1,nT)
#endif
      AeqB = EQ(Coori(1,1),Coori(1,2))
      CeqD = EQ(Coori(1,3),Coori(1,4))
      PrintB10=.False.
      PrintB01=.False.
      PrintB00=.False.
*
      nabMax = la+lb
      ncdMax = lc+ld
      h12 = Half
         If (nabMax.ge.2 .and. ncdMax.ge.2) Then
      Do iT = 1, nT
         Do iRys = 1, nRys
            tmp=h12 * U2(iRys,iT)
            B00(iRys,iT) = tmp
            B10(iRys,iT) = ( h12 - tmp * Eta(iT))*ZInv(iT)
            B01(iRys,iT) = ( h12 - tmp * Zeta(iT))*EInv(iT)
         EndDo
      EndDo
      PrintB10=.True.
      PrintB01=.True.
      PrintB00=.True.
         Else If (ncdMax.eq.0 .and. nabMax.ge.2) Then
      Do iT = 1, nT
         Do iRys = 1, nRys
            B10(iRys,iT) = ( h12 - h12 * U2(iRys,iT) * Eta(iT))*ZInv(iT)
         EndDo
      EndDo
      PrintB10=.True.
         Else If (nabMax.eq.0 .and. ncdMax.ge.2) Then
      Do iT = 1, nT
         Do iRys = 1, nRys
            B01(iRys,iT) =( h12 - h12 * U2(iRys,iT) * Zeta(iT))*EInv(iT)
         EndDo
      EndDo
      PrintB01=.True.
         Else If (ncdMax.eq.1 .and. nabMax.ge.2) Then
      Do iT = 1, nT
         Do iRys = 1, nRys
            tmp=h12 * U2(iRys,iT)
            B00(iRys,iT) = tmp
            B10(iRys,iT) = ( h12 - tmp * Eta(iT))*ZInv(iT)
         EndDo
      EndDo
      PrintB10=.True.
      PrintB00=.True.
         Else If (nabMax.eq.1 .and. ncdMax.ge.2) Then
      Do iT = 1, nT
         Do iRys = 1, nRys
            tmp=h12 * U2(iRys,iT)
            B00(iRys,iT) = tmp
            B01(iRys,iT) = ( h12 - tmp * Zeta(iT))*EInv(iT)
         EndDo
      EndDo
      PrintB01=.True.
      PrintB00=.True.
         Else  If (nabMax.eq.1 .and. ncdMax.eq.1) Then
      Do iT = 1, nT
         Do iRys = 1, nRys
            B00(iRys,iT) = h12*U2(iRys,iT)
         EndDo
      EndDo
      PrintB00=.True.
         End If
*
      If (nabMax.ne.0 .and. ncdMax.ne.0) Then
         If (.Not.AeqB .and. .Not.CeqD) Then
          Do iCar = 1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  tmp=U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
                  PAQP(iRys,iT,iCar) =
     &                  P(iT,iCar) - CoorAC(iCar,1) + Eta(iT)*tmp
                  QCPQ(iRys,iT,iCar) =
     &                  Q(iT,iCar) - CoorAC(iCar,2) - Zeta(iT)*tmp
               EndDo
            EndDo
          EndDo
         Else If (AeqB .and. .Not.CeqD) Then
          Do iCar=1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  tmp=U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
                  PAQP(iRys,iT,iCar) = Eta(iT)*tmp
                  QCPQ(iRys,iT,iCar) =
     &                Q(iT,iCar) - CoorAC(iCar,2) - Zeta(iT)*tmp
               EndDo
            EndDo
          EndDo
         Else If (.Not.AeqB .and. CeqD) Then
          Do iCar=1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  tmp= U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
                  PAQP(iRys,iT,iCar) =
     &                   P(iT,iCar) - CoorAC(iCar,1) + Eta(iT)*tmp
                  QCPQ(iRys,iT,iCar) = - Zeta(iT)*tmp
               EndDo
            EndDo
          EndDo
         Else
          Do iCar=1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  tmp=U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
                   PAQP(iRys,iT,iCar) = Eta(iT)*tmp
                   QCPQ(iRys,iT,iCar) = - Zeta(iT)*tmp
               EndDo
            EndDo
          EndDo
         End If
      Else If (nabMax.ne.0) Then
         If (.Not.AeqB) Then
          Do iCar=1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) =
     &                P(iT,iCar) - CoorAC(iCar,1) + Eta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               EndDo
            EndDo
          EndDo
         Else
          Do iCar=1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  PAQP(iRys,iT,iCar) = Eta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               EndDo
            EndDo
          EndDO
         End If
      Else If (ncdMax.ne.0) Then
         If (.Not.CeqD) Then
          Do iCar=1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  QCPQ(iRys,iT,iCar) =
     &                Q(iT,iCar) - CoorAC(iCar,2) - Zeta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               EndDo
            EndDo
          EndDO
         Else
          Do iCar=1,3
             Do iT = 1, nT
                Do iRys = 1, nRys
                  QCPQ(iRys,iT,iCar) = - Zeta(iT)
     &                * (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
               EndDo
            EndDo
          EndDo
         End If
      End If
#ifdef _DEBUGPRINT_
      If (la+lb.gt.0) Then
         Write (Label,'(A)') ' PAQP(x)'
*        Call RecPrt(Label,' ',PAQP(1,1,1),nRys,nT)
         Write (Label,'(A)') ' PAQP(y)'
*        Call RecPrt(Label,' ',PAQP(1,1,2),nRys,nT)
         Write (Label,'(A)') ' PAQP(z)'
*        Call RecPrt(Label,' ',PAQP(1,1,3),nRys,nT)
      End If
      If (lc+ld.gt.0) Then
         Write (Label,'(A)') ' QCPQ(x)'
*        Call RecPrt(Label,' ',QCPQ(1,1,1),nRys,nT)
         Write (Label,'(A)') ' QCPQ(y)'
*        Call RecPrt(Label,' ',QCPQ(1,1,2),nRys,nT)
         Write (Label,'(A)') ' QCPQ(z)'
*        Call RecPrt(Label,' ',QCPQ(1,1,3),nRys,nT)
      End If
      If (PrintB10) Then
         Write (Label,'(A)') ' B10'
         Call RecPrt(Label,' ',B10(1,1),nRys,nT)
      End If
      If (PrintB00) Then
         Write (Label,'(A)') ' B00'
         Call RecPrt(Label,' ',B00(1,1),nRys,nT)
      End If
      If (PrintB01) Then
         Write (Label,'(A)') ' B01'
         Call RecPrt(Label,' ',B01(1,1),nRys,nT)
      End If
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iDum1)
         Call Unused_integer(iDum2)
         Call Unused_integer(lac)
      End If
      End
