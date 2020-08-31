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
* Copyright (C) 1994, Anders Bernhardsson                              *
*               1994, Roland Lindh                                     *
************************************************************************
      SubRoutine CmbnS2a(Rnxyz,nZeta,la,lb,rKappa,Final,Alpha,
     &                  IfHss,ld)
************************************************************************
*                                                                      *
* Object: compute the 2nd derivative  of the overlap matrix.           *
*                                                                      *
* Called from: OvrHss                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              DDot_   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       rKappa(nZeta),
     &       Rnxyz(nZeta,3,0:la+ld,0:lb), Alpha(nZeta)
      Logical IfHss(4,3,4,3)
      Integer ia(3),ib(3)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
*     Index in the triang. local hessian
*
       iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
      Do 10 iax = 0, la
         ia(1)=iax
         iyaMax=la-ia(1)
      Do 11 ibx = 0, lb
         ib(1)=ibx
         iybMax=lb-ib(1)
         Do 20 iay = 0, iyaMax
            ia(2)=iay
            ia(3) = la-ia(2)-ia(1)
            ipa= Ind(la,ia(1),ia(3))
         Do 21 iby = 0, iybMax
            ib(2)=iby
          ib(3) = lb-ib(2)-ib(1)
          ipb= Ind(lb,ib(1),ib(3))
*
*
*           Combine overlap integrals
*
*           Integrals like dI/dx1dx1
*
          Do 5 iCoor=1,3
            jCoor=Mod(iCoor,3)+1
            kCoor=Mod(jCoor,3)+1
            If (IfHss(1,iCoor,1,iCoor)) Then
                  Do 30 iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,itri(iCoor,iCoor))=rKappa(iZeta)*
     &                   ((Two*Alpha(iZeta))**2 *
     &                   Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor))*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))-
     &                   Two *  Alpha(iZeta) *
     &                   Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                     If (ia(iCoor).gt.0) Then
                        Final(iZeta,ipa,ipb,itri(iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,itri(iCoor,iCoor))
     &                      - rKappa(iZeta)*
     &                      (Four *  Alpha(iZeta)* Dble(ia(iCoor)) *
     &                      Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))*
     &                      Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                      Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                     End If
                     If (ia(iCoor).gt.1) Then
                        Final(iZeta,ipa,ipb,itri(iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,itri(iCoor,iCoor))
     &                      + rKappa(iZeta)*
     &                      (Dble(ia(iCoor)*(ia(iCoor)-1))*
     &                      Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))*
     &                      Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor))*
     &                      Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)))
                     End If
 30               Continue
            End If
 5       Continue
*

*
*           Integrals like dI/dxdz
*
            Do 56 iCoor=2,3
            Do 52 jCoor=1,iCoor-1
                 If (IfHss(1,iCoor,1,jCoor)) Then
                  Do 51 kCoor=1,3
                    Do 50 iZeta = 1, nZeta
                     If (kCoor.eq.1) Then
                        Final(iZeta,ipa,ipb,
     &                       itri(iCoor,jCoor))= rKappa(iZeta)
                     End If
                     If ((kCoor.eq.iCoor).or.(kCoor.eq.jCoor)) Then
                        rIc=Two*Alpha(iZeta)*
     &                      Rnxyz(iZeta,kCoor,ia(kCoor)+1,ib(kCoor))
*
                        If (ia(kCoor).gt.0)
     &                    rIc=rIc-Dble(ia(kCoor))*
     &                       Rnxyz(iZeta,kCoor,ia(kCoor)-1,ib(kCoor))
*
                        Final(iZeta,ipa,ipb,
     &                        itri(iCoor,jCoor))=
     &                  Final(iZeta,ipa,ipb,
     &                        itri(iCoor,jCoor))*
     &                                       rIc
                     Else
                        Final(iZeta,ipa,ipb,
     &                        itri(iCoor,jCoor))=
     &                  Final(iZeta,ipa,ipb,
     &                        itri(iCoor,jCoor))*
     &                                Rnxyz(iZeta,kCoor,ia(kCoor),
     &                                ib(kCoor))
                     End If
 50               Continue
 51             Continue
                End If
 52         Continue
 56         Continue
 21      Continue
 20      Continue
 11   Continue
 10   Continue
      Return
      End
