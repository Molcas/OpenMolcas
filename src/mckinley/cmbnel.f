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
* Copyright (C) 1991,1992,1995, Roland Lindh                           *
*               1997, Anders Bernhardsson                              *
************************************************************************
      SubRoutine CmbnEl(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,Final,nComp,
     &                   Fact,Temp,Alpha,Beta,
     &                   iStb,jStb,nOp,ifgrad,kcar)
************************************************************************
*                                                                      *
* Object: to compute gradient integrals for SC Reaction Fields         *
*                                                                      *
* Called from: RFGrd                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             Modified for reaction field calculations July '92        *
*             Modified for gradient calculations May '95               *
*             Modified for trans. prob.   calculations Oct '97         *
*             by Anders Bernhardsson                                   *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
      Integer nOp(2)
      Logical Ifgrad(3,2)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,2),
     &       Zeta(nZeta), rKappa(nZeta), Fact(nZeta), Temp(nZeta),
     &       Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr),
     &       Alpha(nZeta), Beta(nZeta)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
*
*
      tTwo=Two
*
      Do 130 iZeta = 1, nZeta
         Fact(iZeta) = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two)
130   Continue
*
*---- Loop over angular components of the basis set
*
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 10 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
         Do 20 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
*
*           Combine multipole moment integrals
*
               If (ifgrad(1,1)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (ixa.gt.0) Then
                        xa = -ixa
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   (tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa+1,ixb,ix) +
     &                                xa*Rnxyz(iZeta,1,ixa-1,ixb,ix)) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa+1,ixb,ix) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)-1
                        If (iComp.eq.kcar) Then
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,1) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                        End If
                     End Do
                  End Do
               End Do
               End If
               If (ifgrad(1,2)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (ixb.gt.0) Then
                        xb = -ixb
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   (tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb+1,ix) +
     &                                xb*Rnxyz(iZeta,1,ixa,ixb-1,ix)) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb+1,ix) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)-1
                        If (iComp.eq.kcar) Then
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,2) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                        End If
                     End Do
                  End Do
               End Do
               End If
               If (ifgrad(2,1)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (iya.gt.0) Then
                        ya = -iya
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   (tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,2,iya+1,iyb,iy) +
     &                                ya*Rnxyz(iZeta,2,iya-1,iyb,iy))
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,2,iya+1,iyb,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)-1
                        If (iComp.eq.kcar) Then
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,1) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                        End If
                     End Do
                  End Do
               End Do
               end if
               If (ifgrad(2,2)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (iyb.gt.0) Then
                        yb = -iyb
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   (tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,2,iya,iyb+1,iy) +
     &                                yb*Rnxyz(iZeta,2,iya,iyb-1,iy))
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,2,iya,iyb+1,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)-1
                        If (iComp.eq.kcar) Then
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,2) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                        End If
                     End Do
                  End Do
               End Do
               end if
               If (ifgrad(3,1)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     Do iZeta = 1, nZeta
                        Temp(iZeta) = Fact(iZeta) *
     &                                Rnxyz(iZeta,1,ixa,ixb,ix)*
     &                                Rnxyz(iZeta,2,iya,iyb,iy)
                     End Do
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)-1
                        If (iComp.eq.kcar) Then
                        If (iza.gt.0) Then
                           za = -iza
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,1)= Temp(iZeta)*
     &                                (tTwo*Alpha(iZeta) *
     &                                Rnxyz(iZeta,3,iza+1,izb,iz) +
     &                             za*Rnxyz(iZeta,3,iza-1,izb,iz))
                           End Do
                        Else
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,1)= Temp(iZeta)*
     &                                tTwo*Alpha(iZeta) *
     &                                Rnxyz(iZeta,3,iza+1,izb,iz)
                           End Do
                        End If
                        End If
                     End Do
                  End Do
               End Do
               end if
               If (ifgrad(3,2)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     Do iZeta = 1, nZeta
                        Temp(iZeta) = Fact(iZeta) *
     &                                Rnxyz(iZeta,1,ixa,ixb,ix)*
     &                                Rnxyz(iZeta,2,iya,iyb,iy)
                     End Do
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)-1
                        If (iComp.eq.kcar) Then
                        If (izb.gt.0) Then
                           zb = -izb
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,2)= Temp(iZeta)*
     &                                (tTwo* Beta(iZeta) *
     &                                Rnxyz(iZeta,3,iza,izb+1,iz) +
     &                             zb*Rnxyz(iZeta,3,iza,izb-1,iz))
                           End Do
                        Else
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,2)= Temp(iZeta)*
     &                                tTwo* Beta(iZeta) *
     &                                Rnxyz(iZeta,3,iza,izb+1,iz)
                           End Do
                        End If
                        End If
                     End Do
                  End Do
               End Do
              End If
*
 20      Continue
 10   Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nComp)
         Call Unused_integer(iStb)
         Call Unused_integer(jStb)
         Call Unused_integer_array(nOp)
      End If
      End
