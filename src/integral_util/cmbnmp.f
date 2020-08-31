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
      SubRoutine CmbnMP(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,Final,nComp)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: MltInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), rKappa(nZeta),
     &       Rnxyz(nZeta,3,0:la,0:lb,0:lr)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout = 134
      iPrint = nPrint(iRout)
*     Call qEnter('CmbnMP')
*     Call GetMem(' Enter CmbnMP','LIST','REAL',iDum,iDum)
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
*           If (iPrint.ge.99) Then
*              Write (*,*) ixa,iya,iza,ixb,iyb,izb
*              Write (*,*) ipa,ipb
*           End If
*
*           Combine multipole moment integrals
*
            iComp = 0
            Do 41 ix = lr, 0, -1
               Do 42 iy = lr-ix, 0, -1
                  iz = lr-ix-iy
                  iComp=iComp+1
*                 Write (*,*) ix, iy, iz, iComp
                  Do 30 iZeta = 1, nZeta
                     Fact = rKappa(iZeta) * 1/Sqrt(Zeta(iZeta)**3)
*                    Fact = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two)
                     Final(iZeta,ipa,ipb,iComp) = Fact *
     &                       Rnxyz(iZeta,1,ixa,ixb,ix)*
     &                       Rnxyz(iZeta,2,iya,iyb,iy)*
     &                       Rnxyz(iZeta,3,iza,izb,iz)
 30               Continue
 42            Continue
 41         Continue
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
*
*     Call GetMem(' Exit CmbnMP','LIST','REAL',iDum,iDum)
*     Call qExit('CmbnMP')
      Return
      End
