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
      SubRoutine Util4(nZeta,Final,la,lb,Elalbp,Elalb,Bcoor,Dcoor)
************************************************************************
*                                                                      *
* Object: to assemble the diamagnetic shielding integrals from         *
*         electric field integrals.                                    *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             February '91                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8  Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,9),
     &        Elalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3),
     &        Elalb(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3),
     &        Bcoor(3), Dcoor(3), BD(3)
      Character*80 Label
*
*     Statement function for cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      nElem(ix) = (ix+1)*(ix+2)/2
*
      iRout = 231
      iPrint = nPrint(iRout)
*
      BD(1) = Bcoor(1) - Dcoor(1)
      BD(2) = Bcoor(2) - Dcoor(2)
      BD(3) = Bcoor(3) - Dcoor(3)
      Fact = -1.D6 * One2C2
      If (iPrint.ge.99) Then
          Write (6,*) ' In Util4 la,lb=',la,lb
          Do 100 ia = 1, nElem(la)
             Do 200 ib = 1, nElem(lb+1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbp(',ia,',',ib,',x)'
                Call RecPrt(Label,' ',Elalbp(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbp(',ia,',',ib,',y)'
                Call RecPrt(Label,' ',Elalbp(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbp(',ia,',',ib,',z)'
                Call RecPrt(Label,' ',Elalbp(1,ia,ib,3),nZeta,1)
 200         Continue
 100      Continue
          Do 101 ia = 1, nElem(la)
             Do 201 ib = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalb(',ia,',',ib,',x)'
                Call RecPrt(Label,' ',Elalb(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalb(',ia,',',ib,',y)'
                Call RecPrt(Label,' ',Elalb(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalb(',ia,',',ib,',z)'
                Call RecPrt(Label,' ',Elalb(1,ia,ib,3),nZeta,1)
 201         Continue
 101      Continue
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
               xCxD = Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),1) + BD(1)*
     &              Elalb(iZeta,ipa,ipb,1)
               yCxD = Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2) + BD(1)*
     &              Elalb(iZeta,ipa,ipb,2)
               zCxD = Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3) + BD(1)*
     &              Elalb(iZeta,ipa,ipb,3)
               xCyD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),1) + BD(2)*
     &              Elalb(iZeta,ipa,ipb,1)
               yCyD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),2) + BD(2)*
     &              Elalb(iZeta,ipa,ipb,2)
               zCyD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),3) + BD(2)*
     &              Elalb(iZeta,ipa,ipb,3)
               xCzD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1) + BD(3)*
     &              Elalb(iZeta,ipa,ipb,1)
               yCzD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2) + BD(3)*
     &              Elalb(iZeta,ipa,ipb,2)
               zCzD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3) + BD(3)*
     &              Elalb(iZeta,ipa,ipb,3)
               Final(iZeta,ipa,ipb,1) = Fact * (    +yCyD+zCzD)
               Final(iZeta,ipa,ipb,2) = Fact * (    -yCxD     )
               Final(iZeta,ipa,ipb,3) = Fact * (    -zCxD     )
               Final(iZeta,ipa,ipb,4) = Fact * (    -xCyD     )
               Final(iZeta,ipa,ipb,5) = Fact * (xCxD     +zCzD)
               Final(iZeta,ipa,ipb,6) = Fact * (    -zCyD     )
               Final(iZeta,ipa,ipb,7) = Fact * (    -xCzD     )
               Final(iZeta,ipa,ipb,8) = Fact * (    -yCzD     )
               Final(iZeta,ipa,ipb,9) = Fact * (xCxD+yCyD     )
 30         Continue
*
 21      Continue
 20   Continue
*
 11      Continue
 10   Continue
*
      If (iPrint.ge.49) Then
          Do 300 iComp = 1, 9
             Write (Label,'(A,I2,A)') ' Final (',iComp,') '
             Call RecPrt(Label,' ',Final(1,1,1,iComp),nZeta,
     &                   nElem(la)*nELem(lb))
 300      Continue
      End If
*
      Return
      End
