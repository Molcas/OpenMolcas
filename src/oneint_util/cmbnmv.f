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
      SubRoutine CmbnMV(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,Final,nComp,
     &                  rV2Int,rV4Int)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: MVeInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN, February '91                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,nComp,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2),
     &       Zeta(nZeta), rKappa(nZeta),
     &       Rnxyz(nZeta,3,0:la+2,0:lb+2,0:lr),
     &       rV2Int(nZeta,3,0:la,0:lb,2),
     &       rV4Int(nZeta,3,0:la,0:lb)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout = 191
      iPrint = nPrint(iRout)
*     Call GetMem(' Enter CmbnMV','LIST','REAL',iDum,iDum)
*
      Const = - One2C2/Four
      iComp = 1
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
*           Combine integrals
*
            Do 30 iZeta = 1, nZeta
               Fact = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two) * Const
               x2x2 =   rV4Int(iZeta,1,ixa,ixb)*
     &                  Rnxyz(iZeta,2,iya,iyb,0)*
     &                  Rnxyz(iZeta,3,iza,izb,0)
               x2y2 =   rV2Int(iZeta,1,ixa,ixb,1)*
     &                  rV2Int(iZeta,2,iya,iyb,2)*
     &                  Rnxyz(iZeta,3,iza,izb,0)
               x2z2 =   rV2Int(iZeta,1,ixa,ixb,1)*
     &                  Rnxyz(iZeta,2,iya,iyb,0)*
     &                  rV2Int(iZeta,3,iza,izb,2)
               y2x2 =   rV2Int(iZeta,1,ixa,ixb,2)*
     &                  rV2Int(iZeta,2,iya,iyb,1)*
     &                  Rnxyz(iZeta,3,iza,izb,0)
               y2y2 =   Rnxyz(iZeta,1,ixa,ixb,0)*
     &                  rV4Int(iZeta,2,iya,iyb)*
     &                  Rnxyz(iZeta,3,iza,izb,0)
               y2z2 =   Rnxyz(iZeta,1,ixa,ixb,0)*
     &                  rV2Int(iZeta,2,iya,iyb,1)*
     &                  rV2Int(iZeta,3,iza,izb,2)
               z2x2 =   rV2Int(iZeta,1,ixa,ixb,2)*
     &                  Rnxyz(iZeta,2,iya,iyb,0)*
     &                  rV2Int(iZeta,3,iza,izb,1)
               z2y2 =   Rnxyz(iZeta,1,ixa,ixb,0)*
     &                  rV2Int(iZeta,2,iya,iyb,2)*
     &                  rV2Int(iZeta,3,iza,izb,1)
               z2z2 =   Rnxyz(iZeta,1,ixa,ixb,0)*
     &                  Rnxyz(iZeta,2,iya,iyb,0)*
     &                  rV4Int(iZeta,3,iza,izb)
*
               rMVel= x2x2+x2y2+x2z2+y2x2+y2y2+y2z2+z2x2+z2y2+z2z2
*
               Final(iZeta,iComp,ipa,ipb) = Fact * rMVel
 30         Continue
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
*
*     Call GetMem(' Exit CmbnMV','LIST','REAL',iDum,iDum)
      Return
      End
