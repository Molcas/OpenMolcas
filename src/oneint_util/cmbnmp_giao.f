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
* Copyright (C) 1991,2000, Roland Lindh                                *
************************************************************************
      SubRoutine CmbnMP_GIAO(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,Final,
     &                       nComp,nB,RAB,C)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*                                                                      *
*             Modified to GIAO 1st derivatives by R. Lindh in Tokyo,   *
*              Japan, January 2000.                                    *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp,nB),
     &       Zeta(nZeta), rKappa(nZeta),
     &       Rnxyz(nZeta,3,0:la,0:lb,0:lr+1), RAB(3), C(3)
      Integer ix_(3,2)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
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
*           Combine multipole moment integrals
*
            Do iBx = 1, 3
               iBy = iBx + 1
               If (iBy.gt.3) iBy=iBy-3
               iBz = iBy + 1
               If (iBz.gt.3) iBz=iBz-3
               Call ICopy(6,[0],0,ix_,1)
               ix_(iBz,1)=1
               ix_(iBy,2)=1
            iComp = 0
            Do ix = lr, 0, -1
               Do iy = lr-ix, 0, -1
                  iz = lr-ix-iy
                  iComp=iComp+1
                  Do iZeta = 1, nZeta
                     Fact = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two)
                     temp =  Rnxyz(iZeta,1,ixa,ixb,ix)*
     &                       Rnxyz(iZeta,2,iya,iyb,iy)*
     &                       Rnxyz(iZeta,3,iza,izb,iz)
                     tempz=  Rnxyz(iZeta,1,ixa,ixb,ix+ix_(1,1))*
     &                       Rnxyz(iZeta,2,iya,iyb,iy+ix_(2,1))*
     &                       Rnxyz(iZeta,3,iza,izb,iz+ix_(3,1))
                     tempy=  Rnxyz(iZeta,1,ixa,ixb,ix+ix_(1,2))*
     &                       Rnxyz(iZeta,2,iya,iyb,iy+ix_(2,2))*
     &                       Rnxyz(iZeta,3,iza,izb,iz+ix_(3,2))
*
*------------------- The term has only an imaginary component
*
                     Final(iZeta,ipa,ipb,iComp,iBx) = Half * Fact * (
     &                       RAB(iBy)*(Tempz+C(iBz)*Temp)
     &                    -  RAB(iBz)*(Tempy+C(iBy)*Temp)
     &                                                      )
                  End Do
               End Do
            End Do
            End Do  ! iB
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
*
      Return
      End
