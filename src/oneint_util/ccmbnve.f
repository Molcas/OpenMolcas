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
      SubRoutine CCmbnVe(Rnxyz,nZeta,la,lb,Zeta,rKappa,Final,nComp,
     &                  Vxyz,KVector)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January  91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), rKappa(nZeta), rTemp, Fact, KVector(3)
      Complex*16 Rnxyz(nZeta,3,0:la+1,0:lb+1),
     &          Vxyz(nZeta,3,0:la,0:lb,2), Temp1, Temp2
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout = 161
      iPrint = nPrint(iRout)
*
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 11 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
         Do 21 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
*
*           Combine integrals
*
            Do 30 iZeta = 1, nZeta
*
*              Put in the correct prefactors
*
               rTemp=KVector(1)**2 + kVector(2)**2 + kVector(3)**2
               rTemp=rTemp/(Four*Zeta(iZeta))
               Fact = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two) *
     &               Exp(-rTemp)
               Temp1= Fact *
     &               Vxyz(iZeta,1,ixa,ixb,1)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb)
               Temp2= Fact *
     &               Vxyz(iZeta,1,ixa,ixb,2)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb)
               Final(iZeta,ipa,ipb,1) = DBLE((Temp1+Temp2)*Half)
               Final(iZeta,ipa,ipb,4) = DBLE((Temp1-Temp2)*Half)
               Final(iZeta,ipa,ipb,7) = DIMAG((Temp1+Temp2)*Half)
               Final(iZeta,ipa,ipb,10)= DIMAG((Temp1-Temp2)*Half)
               Temp1= Fact *
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &               Vxyz(iZeta,2,iya,iyb,1)*
     &               Rnxyz(iZeta,3,iza,izb)
               Temp2= Fact *
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &               Vxyz(iZeta,2,iya,iyb,2)*
     &               Rnxyz(iZeta,3,iza,izb)
               Final(iZeta,ipa,ipb,2) = DBLE((Temp1+Temp2)*Half)
               Final(iZeta,ipa,ipb,5) = DBLE((Temp1-Temp2)*Half)
               Final(iZeta,ipa,ipb,8) = DIMAG((Temp1+Temp2)*Half)
               Final(iZeta,ipa,ipb,11)= DIMAG((Temp1-Temp2)*Half)
               Temp1= Fact *
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Vxyz(iZeta,3,iza,izb,1)
               Temp2= Fact *
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Vxyz(iZeta,3,iza,izb,2)
               Final(iZeta,ipa,ipb,3) = DBLE((Temp1+Temp2)*Half)
               Final(iZeta,ipa,ipb,6 )= DBLE((Temp1-Temp2)*Half)
               Final(iZeta,ipa,ipb,9 )= DIMAG((Temp1+Temp2)*Half)
               Final(iZeta,ipa,ipb,12)= DIMAG((Temp1-Temp2)*Half)
 30         Continue
            If (iPrint.ge.99) Then
               Write (6,*) '(',ixa,iya,iza,ixb,iyb,izb,')'
               Write (6,*) 'x-component'
               Write (6,*) Final(1,ipa,ipb,1)
               Write (6,*) Final(1,ipa,ipb,4)
               Write (6,*) Final(1,ipa,ipb,7)
               Write (6,*) Final(1,ipa,ipb,10)
               Write (6,*) 'y-component'
               Write (6,*) Final(1,ipa,ipb,2)
               Write (6,*) Final(1,ipa,ipb,5)
               Write (6,*) Final(1,ipa,ipb,8)
               Write (6,*) Final(1,ipa,ipb,11)
               Write (6,*) 'z-component'
               Write (6,*) Final(1,ipa,ipb,3)
               Write (6,*) Final(1,ipa,ipb,6)
               Write (6,*) Final(1,ipa,ipb,9)
               Write (6,*) Final(1,ipa,ipb,12)
            End If
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
*
      Return
      End
