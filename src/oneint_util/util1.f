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
      SubRoutine Util1(Alpha,Beta,nZeta,Final,la,lb,
     *                 Slaplb,Slamlb,Slalbp,Slalbm)
************************************************************************
*                                                                      *
* Object: to assemble the electric field integrals from                *
*         derivative integrals of the electric potential.              *
*                                                                      *
* Called from: EFInt                                                   *
*                                                                      *
* Calling    : qEnter                                                  *
*              qExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             February '91                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8  Final(nZeta,3,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2),
     *        Slaplb(nZeta,(la+2)*(la+3)/2,(lb+1)*(lb+2)/2),
     *        Slamlb(nZeta,(la)*(la+1)/2,(lb+1)*(lb+2)/2),
     *        Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2),
     *        Slalbm(nZeta,(la+1)*(la+2)/2,(lb)*(lb+1)/2),
     *        Alpha(nZeta), Beta(nZeta)
      Character*80 Label
*
*     Statement function for cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      nElem(ix) = (ix+1)*(ix+2)/2
*
      iRout = 203
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
          Write (6,*) ' In Util1 la,lb=',la,lb
          Call RecPrt('Alpha',' ',Alpha,nZeta,1)
          Call RecPrt('Beta',' ',Beta,nZeta,1)
          Do 200 ib = 1, nElem(lb)
             Write (Label,'(A,I2,A)') ' Slaplb(la,',ib,')'
             Call RecPrt(Label,' ',Slaplb(1,1,ib),nZeta,nElem(la+1))
200       Continue
          If (la.gt.0) Then
             Do 201 ib = 1, nElem(lb)
                Write (Label,'(A,I2,A)') ' Slamlb(la,',ib,')'
                Call RecPrt(Label,' ',Slamlb(1,1,ib),nZeta,nElem(la-1))
201          Continue
          End If
          Do 300 ib = 1, nElem(lb+1)
             Write (Label,'(A,I2,A)') ' Slalbp(la,',ib,')'
             Call RecPrt(Label,' ',Slalbp(1,1,ib),nZeta,nElem(la))
300       Continue
          If (lb.gt.0) Then
             Do 301 ib = 1, nElem(lb-1)
                Write (Label,'(A,I2,A)') ' Slalbm(la,',ib,')'
                Call RecPrt(Label,' ',Slalbm(1,1,ib),nZeta,nElem(la))
301          Continue
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
            If (ixa.eq.0 .and. ixb.eq.0) Then
*
            Do 33 iZeta = 1, nZeta
               Final(iZeta,1,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
33          Continue
*
            Else If (ixa.eq.0) Then
*
            Do 32 iZeta = 1, nZeta
               Final(iZeta,1,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
     *          -Dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
32          Continue
*
            Else If (ixb.eq.0) Then
*
            Do 31 iZeta = 1, nZeta
               Final(iZeta,1,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)
     *          -Dble(ixa)*Slamlb(iZeta,Ind(la-1,ixa-1,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
31          Continue
*
            Else
*
            Do 30 iZeta = 1, nZeta
               Final(iZeta,1,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)
     *          -Dble(ixa)*Slamlb(iZeta,Ind(la-1,ixa-1,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
     *          -Dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
30          Continue
*
            End If
*
            If (iya.eq.0 .and. iyb.eq.0) Then
*
            Do 63 iZeta = 1, nZeta
               Final(iZeta,2,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
63          Continue
*
            Else If (iya.eq.0) Then
*
            Do 62 iZeta = 1, nZeta
               Final(iZeta,2,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
     *          -Dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
62          Continue
*
            Else If (iyb.eq.0) Then
*
            Do 61 iZeta = 1, nZeta
               Final(iZeta,2,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)
     *          -Dble(iya)*Slamlb(iZeta,Ind(la-1,ixa,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
61          Continue
*
            Else
*
            Do 60 iZeta = 1, nZeta
               Final(iZeta,2,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)
     *          -Dble(iya)*Slamlb(iZeta,Ind(la-1,ixa,iza),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
     *          -Dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
60          Continue
*
            End If
*
            If (iza.eq.0 .and. izb.eq.0) Then
*
            Do 93 iZeta = 1, nZeta
               Final(iZeta,3,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
93          Continue
*
            Else If (iza.eq.0) Then
*
            Do 92 iZeta = 1, nZeta
               Final(iZeta,3,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
     *          -Dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
92          Continue
*
            Else If (izb.eq.0) Then
*
            Do 91 iZeta = 1, nZeta
               Final(iZeta,3,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)
     *          -Dble(iza)*Slamlb(iZeta,Ind(la-1,ixa,iza-1),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
91          Continue
*
            Else
*
            Do 90 iZeta = 1, nZeta
               Final(iZeta,3,ipa,ipb) =
     *           Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)
     *          -Dble(iza)*Slamlb(iZeta,Ind(la-1,ixa,iza-1),ipb)
     *          +Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
     *          -Dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
90          Continue
*
            End If
*
21       Continue
20    Continue
*
11       Continue
10    Continue
*
      If (iPrint.ge.49) Then
          Write (6,*) ' In Util1 la,lb=',la,lb
          Do 400 iElem = 1, nElem(la)
             Do 410 jElem = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)')
     *                ' Final (',iElem,',',jElem,') '
                Call RecPrt(Label,' ',Final(1,1,iElem,jElem),nZeta,3)
410          Continue
400       Continue
      End If
*
      Return
      End
