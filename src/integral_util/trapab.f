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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine TraPAB(nZeta,la,lb,AB,GInt,jSum,rKappa,Fac1,Fac2,
     &                  Fac3,Fac4,Fac5,A,B,P)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
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
      Real*8 AB(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2),GInt(nZeta,jSum),
     &       rKappa(nZeta), Fac1(nZeta), Fac2(nZeta), Fac3(nZeta),
     &       Fac4(nZeta), Fac5(nZeta), A(3), B(3), P(nZeta,3)
*
*-----Statement function
*
      iad(ix,iy,iz) = (iy+iz) * (iy+iz+1) / 2 + iz + 1
      iOff(ix,iy,iz) = (ix+iy+iz)*(ix+iy+iz+1)*(ix+iy+iz+2)/6
      nElem(i) = (i+1)*(i+2)/2
*
      iRout = 239
      iPrint = nPrint(iRout)
      Call QEnter('TraPAB')
      If (iPrint.ge.99) Then
         Call RecPrt(' In TraPAB: GInt',' ',GInt,nZeta,jSum)
         Call RecPrt(' In TraPAB: P   ',' ',P   ,nZeta,3)
      End If
*
*-----Initilize
*
      call dcopy_(nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2,Zero,0,AB,1)
*
*-----Remove redundant elements in GInt. This is done in place.
*
      kOff = 4
      Do 101 i = 2, la+lb
*
         Do 110 ix = i, 0, -1
            Do 120 iy = i-ix, 0, -1
               iz = i-ix-iy
               jx = ix
               jy = iy
               jz = iz
*
               lOff = 0
               Do 121 ia = 1, i-1
                  If (jz.ne.0) Then
                     lOff = 3*(2+lOff)
                     jz = jz - 1
                  Else If (jy.ne.0) Then
                     lOff = 3*(1+lOff)
                     jy = jy - 1
                  Else
                     lOff = 3*lOff
                     jx = jx - 1
                  End If
 121           Continue
               If (jz.eq.1) lOff = lOff + 3
               If (jy.eq.1) lOff = lOff + 2
               If (jx.eq.1) lOff = lOff + 1
*
               iTrgt = iOff(ix,iy,iz) + iAd(ix,iy,iz)
*              Write(*,*) ' ix,iy,iz,kOff,lOff,iTrgt=',
*    &                      ix,iy,iz,kOff,lOff,iTrgt
               call dcopy_(nZeta,GInt(1,kOff+lOff),1,GInt(1,iTrgt),1)
*
 120        Continue
 110     Continue
*
         kOff = kOff + 3**i
 101  Continue
      If (iPrint.ge.99) Then
         Call RecPrt(' In TraPAB: GInt(unique)',' ',GInt,nZeta,
     &      (la+lb+1)*(la+lb+2)*(la+lb+3)/6)
      End If
*
*-----Loop over the elements of the basis functions on A and B
*
      Do 10 ixa = la, 0, -1
         iyaMax = la - ixa
         Do 20 iya = iyaMax, 0, -1
            iza = la - ixa - iya
            ipa = iad(ixa,iya,iza)
*           Write (*,*) ' ipa,ixa,iya,iza=',ipa,ixa,iya,iza
*
            Do 40 ixb = lb, 0, -1
               iybMax = lb - ixb
               Do 50 iyb = iybMax, 0, -1
                  izb = lb - ixb - iyb
                  ipb = iad(ixb,iyb,izb)
*                 Write (*,*) ' ipb,ixb,iyb,izb=',ipb,ixb,iyb,izb
*
*-----------------Loop over the elements of functions at P
*
                  Do 11 ixas = 0, ixa
                     Call Binom(ixa,ixas,iAx)
                     Ax = DBLE(iAx)
                     Do 12 iZeta = 1, nZeta
                       If ( (ixa-ixas).eq.0 ) then
                        Fac1(iZeta) = rKappa(iZeta) * Ax
                       Else
                        Fac1(iZeta) = rKappa(iZeta) * Ax *
     &                      (P(iZeta,1)-A(1))**(ixa-ixas)
                       End If
 12                  Continue
                     Do 21 iyas = 0, iya
                        Call Binom(iya,iyas,iAy)
                        Ay = DBLE(iAy)
                        Do 22 iZeta = 1, nZeta
                          If ( (iya-iyas).eq.0 ) then
                           Fac2(iZeta) = Fac1(iZeta) * Ay
                          Else
                           Fac2(iZeta) = Fac1(iZeta) * Ay *
     &                         (P(iZeta,2)-A(2))**(iya-iyas)
                          End If
 22                     Continue
                        Do 31 izas = 0, iza
                           Call Binom(iza,izas,iAz)
                           Az = DBLE(iAz)
                           Do 32 iZeta = 1, nZeta
                             If ( (iza-izas).eq.0 ) then
                              Fac3(iZeta) = Fac2(iZeta) * Az
                             Else
                              Fac3(iZeta) = Fac2(iZeta) * Az *
     &                            (P(iZeta,3)-A(3))**(iza-izas)
                             End If
 32                        Continue
*
                           Do 41 ixbs = 0, ixb
                              Call Binom(ixb,ixbs,iBx)
                              Bx = DBLE(iBx)
                              Do 42 iZeta = 1, nZeta
                                If ( (ixb-ixbs).eq.0 ) then
                                 Fac4(iZeta) = Fac3(iZeta) * Bx
                                Else
                                 Fac4(iZeta) = Fac3(iZeta) * Bx *
     &                               (P(iZeta,1)-B(1))**(ixb-ixbs)
                                End If
 42                           Continue
                              igx = ixas + ixbs
                              Do 51 iybs = 0, iyb
                                 Call Binom(iyb,iybs,iBy)
                                 By = DBLE(iBy)
                                 Do 52 iZeta = 1, nZeta
                                   If ( (iyb-iybs).eq.0 ) then
                                    Fac5(iZeta) = Fac4(iZeta) * By
                                   Else
                                    Fac5(iZeta) = Fac4(iZeta) * By *
     &                                  (P(iZeta,2)-B(2))**(iyb-iybs)
                                   End If
 52                              Continue
                                 igy = iyas + iybs
                                 Do 61 izbs = 0, izb
                                    Call Binom(izb,izbs,iBz)
                                    Bz = DBLE(iBz)
                                    igz = izas + izbs
                                    ipg = iOff(igx,igy,igz) +
     &                                    iAd(igx,igy,igz)
*                 Write (*,*) ' ipg,igx,igy,igz=', ipg,igx,igy,igz
*
                                    Do 100 iZeta = 1, nZeta
                                      If ( (izb-izbs).eq.0 ) then
                                       AB(iZeta,ipa,ipb) =
     &                                     AB(iZeta,ipa,ipb) +
     &                                     Fac5(iZeta) *
     &                                     GInt(iZeta,ipg) * Bz
                                      Else
                                       AB(iZeta,ipa,ipb) =
     &                                     AB(iZeta,ipa,ipb) +
     &                                     Fac5(iZeta) *
     &                                     (P(iZeta,3)-B(3))**(izb-izbs)
     &                                     * GInt(iZeta,ipg) * Bz
                                      End If
 100                                Continue
*
 61                              Continue
 51                           Continue
 41                        Continue
*
 31                     Continue
 21                  Continue
 11               Continue
*
 50            Continue
 40         Continue
*
 20      Continue
 10   Continue
*
      If (iPrint.ge.89) Then
         nab = nElem(la) * nElem(lb)
         Call RecPrt(' In TraPAB: AB',' ', AB, nZeta,nab)
      End If
*
      Call QExit('TraPAB')
      Return
      End
