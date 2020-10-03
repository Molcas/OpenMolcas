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
* Copyright (C) 1994, Bernd Artur Hess                                 *
************************************************************************
      SubRoutine util8(Beta,nZeta,Final,la,lb,Slalbp,Slalbm)
************************************************************************
*                                                                      *
* Object: to assemble the Vp integrals from                            *
*         derivative integrals of the electric potential.              *
*                                                                      *
*     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
*             Chemie, University of Bonn, Germany, August 1994         *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8  Final (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3),
     *        Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2),
     *        Slalbm(nZeta,(la+1)*(la+2)/2, lb   *(lb+1)/2),
     *        Beta(nZeta)
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
*
      If (iPrint.ge.99) Then
          Write (6,*) ' In util8 la,lb=',la,lb
          Call RecPrt('Beta','(5f15.8)',Beta,nZeta,1)
          Do 200 ib = 1, nElem(lb)
             Write (Label,'(A,I2,A)') ' Slalbp(',la,ib,')'
      Call RecPrt(Label,'(5f15.8)',Slalbp(1,1,ib),nZeta,nElem(la+1))
200       Continue
          If (lb.gt.0) Then
             Do 201 ia = 1, nElem(la)
                Write (Label,'(A,I2,A)') ' Slalbm(',la,ib,')'
      Call RecPrt(Label,'(5f15.8)',Slalbm(1,1,ib),nZeta,nElem(lb-1))
201          Continue
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
            IF (ixb.EQ.0) THEN
             Do 33 iZeta = 1, nZeta
                Final(iZeta,ipa,ipb,1) =
     *         Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
33           Continue
            ELSE
             Do 31 iZeta = 1, nZeta
                Final(iZeta,ipa,ipb,1) =
     *         Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
     *           -Dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
31           Continue
            ENDIF
*
            IF (iyb.EQ.0) THEN
*
            Do 63 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,2) =
     *           Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
63          Continue
            ELSE
            Do 62 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,2) =
     *           Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
     *          -Dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
62          Continue
            ENDIF
*
            IF (izb.EQ.0) THEN
*
            Do 93 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,3) =
     *         Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
93          Continue
            ELSE
            Do 91 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,3) =
     *         Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
     *          -Dble(iza)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
91          Continue
            ENDIF
*
21       Continue
20    Continue
*
11       Continue
10    Continue
*
      If (iPrint.ge.49) Then
          Write (6,*) ' In UTIL8 la,lb=',la,lb
          Do 380 iiComp=1,3
             Do 410 jElem = 1, nElem(lb)
          Do 400 iElem = 1, nElem(la)
          Do 390 iiZeta=1,nZeta
                Write (Label,'(A,I2,A,I2,A)')
     *                ' Final (',iElem,',',jElem,') '
*      Call RecPrt(Label,'(5f15.8)',
         write(6,*) iiZeta,iElem,jElem,iiComp,
     *        Final(iiZeta,iElem,jElem,iiComp)
390        Continue
400       Continue
410          Continue
380          Continue
      End If
*
      Return
      End
