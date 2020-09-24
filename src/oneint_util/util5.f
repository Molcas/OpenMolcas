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
*               2019, Ignacio Fdez. Galvan                             *
************************************************************************
      SubRoutine Util5(Beta,nZeta,Final,la,lb,Slalbp,Slalbm)
************************************************************************
*                                                                      *
* Object: to assemble the velocity quadrupole integrals from the       *
*         derivative integrals and dipole integrals.                   *
*                                                                      *
* Called from: QpVInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             February '91                                             *
*     Adapted from util2: Ignacio Fdez. Galvan, July 2019              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8  Final (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &        Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3),
     &        Slalbm(nZeta,(la+1)*(la+2)/2, lb   *(lb+1)/2,3),
     &        Beta(nZeta)
      Character*80 Label
*
*     Statement function for cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      nElem(ix) = (ix+1)*(ix+2)/2
*
      iRout = 211
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
          Write (6,*) ' In Util5 la,lb=',la,lb
          Call RecPrt('Beta',' ',Beta,nZeta,1)
          Do 100 ia = 1, nElem(la)
             Do 200 ib = 1, nElem(lb+1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbp(',ia,',',ib,'x)'
                Call RecPrt(Label,' ',Slalbp(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbp(',ia,',',ib,'y)'
                Call RecPrt(Label,' ',Slalbp(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbp(',ia,',',ib,'z)'
                Call RecPrt(Label,' ',Slalbp(1,ia,ib,3),nZeta,1)
 200         Continue
 100      Continue
          If (lb.gt.0) Then
             Do 101 ia = 1, nElem(la)
                Do 201 ib = 1, nElem(lb-1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Slalbm(',ia,',',ib,'x)'
                   Call RecPrt(Label,' ',Slalbm(1,ia,ib,1),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Slalbm(',ia,',',ib,'y)'
                   Call RecPrt(Label,' ',Slalbm(1,ia,ib,2),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Slalbm(',ia,',',ib,'z)'
                   Call RecPrt(Label,' ',Slalbm(1,ia,ib,3),nZeta,1)
 201            Continue
 101         Continue
          End If
      End If
*
      Do 10 ixa = la, 0, -1
         Do 11 iya = la-ixa, 0, -1
            iza = la-ixa-iya
            ipa = Ind(la,ixa,iza)
*
      Do 20 ixb = lb, 0, -1
         Do 21 iyb = lb-ixb, 0, -1
            izb = lb-ixb-iyb
            ipb = Ind(lb,ixb,izb)
*
            Do 30 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,1) = Two*Beta(iZeta)*
     &            Two*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),1)
               Final(iZeta,ipa,ipb,2) = Two*Beta(iZeta)*
     &            (Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb)  ,1)+
     &             Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2))
               Final(iZeta,ipa,ipb,3) = Two*Beta(iZeta)*
     &            (Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1)+
     &             Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3))
               Final(iZeta,ipa,ipb,4) = Two*Beta(iZeta)*
     &            Two*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb)  ,2)
               Final(iZeta,ipa,ipb,5) = Two*Beta(iZeta)*
     &            (Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2)+
     &             Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb)  ,3))
               Final(iZeta,ipa,ipb,6) = Two*Beta(iZeta)*
     &            Two*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3)
 30         Continue
*
            If (ixb.gt.0) Then
               Do 31 iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,1) = Final(iZeta,ipa,ipb,1)
     &                -Dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),1)
     &                 *Two
                  Final(iZeta,ipa,ipb,2) = Final(iZeta,ipa,ipb,2)
     &                -Dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),2)
                  Final(iZeta,ipa,ipb,3) = Final(iZeta,ipa,ipb,3)
     &                -Dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),3)
 31            Continue
            End If
*
            If (iyb.gt.0) Then
               Do 32 iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,2) = Final(iZeta,ipa,ipb,2)
     &                -Dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb)  ,1)
                  Final(iZeta,ipa,ipb,4) = Final(iZeta,ipa,ipb,4)
     &                -Dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb)  ,2)
     &                 *Two
                  Final(iZeta,ipa,ipb,5) = Final(iZeta,ipa,ipb,5)
     &                -Dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb)  ,3)
 32            Continue
            End If
*
            If (izb.gt.0) Then
               Do 33 iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,3) = Final(iZeta,ipa,ipb,3)
     &                -Dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),1)
                  Final(iZeta,ipa,ipb,5) = Final(iZeta,ipa,ipb,5)
     &                -Dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2)
                  Final(iZeta,ipa,ipb,6) = Final(iZeta,ipa,ipb,6)
     &                -Dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3)
     &                 *Two
 33            Continue
            End If
*
 21      Continue
 20   Continue
*
 11      Continue
 10   Continue
*
      If (iPrint.ge.49) Then
          Write (6,*) ' In Util5 la,lb=',la,lb
          Do 300 iElem = 1, nElem(la)
             Do 310 jElem = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',xx) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',xy) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',xz) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,3),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',yy) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,4),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',yz) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,5),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',zz) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,6),nZeta,1)
 310         Continue
 300      Continue
      End If
*
      Return
      End
