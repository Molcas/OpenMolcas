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
      SubRoutine XCff2D(iDum1,iDum2,nRys,
     &                  Zeta,ZInv,rDum3,rDum4,nT,
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
      Real*8 Zeta(nT), ZInv(nT),
     &       Coori(3,4), CoorAC(3,2),
     &       P(nT,3), Q(nT,3), U2(nRys,nT),
     &       PAQP(nRys,nT,3), QCPQ(nRys,nT,3),
     &       B10(nRys,nT,3),
     &       B00(nRys,nT,3),
     &       B01(nRys,nT,3)
*     Local arrays
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Character*30 Label
#endif
      Logical AeqB, CeqD, EQ, lB10, lB00, lB01
*
      iRout = 14
      iPrint = nPrint(iRout)
*
#ifdef _DEBUGPRINT_
      iPrint=99
      If (iPrint.ge.99) Then
         Call RecPrt(' In XCff2D: Coori',' ',Coori,3,4)
         Call RecPrt(' In XCff2D: P',' ',P,nT,3)
         Call RecPrt(' In XCff2D: Q',' ',Q,nT,3)
      End If
#endif
      AeqB = EQ(Coori(1,1),Coori(1,2))
      CeqD = EQ(Coori(1,3),Coori(1,4))
*
      nabMax=la+lb
      ncdMax=ld+lc
      h12 = Half
*
*---- Compute B10, B00, and B01
*
      lB10=.False.
      lB00=.False.
      lB01=.False.
      If (nabMax.gt.1) Then
         lB10=.True.
         Do 10 iT = 1, nT
            Do 31 iRys = 1, nRys
                  B10(iRys,iT,1) = ( h12 -
     &              h12 * U2(iRys,iT))*ZInv(iT)
 31         Continue
 10      Continue
         call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,2),1)
         call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,3),1)
      End If
      If (lac.ne.0) Then
         lB00=.True.
         call dcopy_(nRys*nT,U2(1,1),1,B00(1,1,1),1)
         call dcopy_(nRys*nT,U2(1,1),1,B00(1,1,2),1)
         call dcopy_(nRys*nT,U2(1,1),1,B00(1,1,3),1)
      End If
      If (ncdMax.gt.1) Then
         lB01=.True.
         Do iT = 1, nT
            Do iRys = 1, nRys
               B01(iRys,iT,1) = Two * Zeta(iT) * U2(iRys,iT)
            End Do
         End Do
         call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,2),1)
         call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,3),1)
      End If
*
      If (nabMax.ne.0 .and. ncdMax.ne.0) Then
         If (.Not.AeqB .and. CeqD) Then
            Do 300 iCar = 1, 3
               Do 310 iT = 1, nT
                  Do 330 iRys = 1, nRys
                     PAQP(iRys,iT,iCar) =
     &                   P(iT,iCar) - CoorAC(iCar,1) +
     &                   (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
                     QCPQ(iRys,iT,iCar) = - Two * Zeta(iT) *
     &                   (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
 330              Continue
 310           Continue
 300        Continue
         Else
            Do 400 iCar = 1, 3
               Do 410 iT = 1, nT
                  Do 430 iRys = 1, nRys
                     PAQP(iRys,iT,iCar) =
     &                   (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
                     QCPQ(iRys,iT,iCar) = - Two * Zeta(iT) *
     &                   (U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar)))
 430              Continue
 410           Continue
 400        Continue
         End If
      Else If (nabMax.ne.0) Then
         If (.Not.AeqB) Then
            Do 101 iCar = 1, 3
               Do 111 iT = 1, nT
                  Do 131 iRys = 1, nRys
                     PAQP(iRys,iT,iCar) =
     &                   P(iT,iCar) - CoorAC(iCar,1) +
     &                   U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
 131              Continue
 111           Continue
 101        Continue
         Else
            Do 201 iCar = 1, 3
               Do 211 iT = 1, nT
                  Do 231 iRys = 1, nRys
                     PAQP(iRys,iT,iCar) =
     &                    U2(iRys,iT) * (Q(iT,iCar)-P(iT,iCar))
 231              Continue
 211           Continue
 201        Continue
         End If
      Else If (ncdMax.ne.0) Then
         Do 202 iCar = 1, 3
            Do 212 iT = 1, nT
               Do 232 iRys = 1, nRys
                  QCPQ(iRys,iT,iCar) = Two * Zeta(iT) *
     &                U2(iRys,iT) * (P(iT,iCar)-Q(iT,iCar))
 232              Continue
 212        Continue
 202     Continue
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
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iDum1)
         Call Unused_integer(iDum2)
         Call Unused_real(rDum3)
         Call Unused_real(rDum4)
      End If
      End
