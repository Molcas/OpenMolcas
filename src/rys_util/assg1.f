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
*               1996, Hans-Joachim Werner                              *
************************************************************************
      SubRoutine Assg1(Temp,PAO,nT,nRys,la,lb,lc,ld,xyz2D0,xyz2D1,
     &                 IfGrad,Index,mVec)
************************************************************************
*                                                                      *
* Object: to assemble the gradients of the ERI's.                      *
*                                                                      *
* Called from: Rysg1                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91; modified by H.-J. Werner, Mai 1996          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "iavec.fh"
      Real*8 PAO(nT,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,
     &              (lc+1)*(lc+2)/2,(ld+1)*(ld+2)/2),
     &       xyz2D0(nRys,nT,0:la+1,0:lb+1,0:lc+1,0:ld+1,3),
     &       xyz2D1(nRys,nT,0:la  ,0:lb  ,0:lc  ,0:ld  ,9),
     &      Temp(9)
      Logical IfGrad(3,4)
      Integer Ind1(3,3), Ind2(3,3), Index(3,4), nVec(3)
*
*     Statement functions
*
      nElem(i) = (i+1)*(i+2)/2
*
      iRout = 248
      iPrint = nPrint(iRout)
*
      call dcopy_(9,Zero,0,Temp,1)
*
      ii = la*(la+1)*(la+2)/6
      jj = lb*(lb+1)*(lb+2)/6
      kk = lc*(lc+1)*(lc+2)/6
      ll = ld*(ld+1)*(ld+2)/6
*
      mVec = 0
      Do i = 1, 3             ! Cartesian directions
         nVec(i) = 0
         Do iCent = 1, 4      ! Centers of integral
            If (IfGrad(i,iCent)) Then
               mVec = mVec + 1
               nVec(i) = nVec(i) + 1
               Ind1(nVec(i),i) = 3*(Index(i,iCent)-1) + i
               Ind2(nVec(i),i) = mVec
            End If
         End Do
      End Do
*
      Do 100 ipd = 1, nElem(ld)
         ixd = ixyz(1,ll+ipd)
         iyd = ixyz(2,ll+ipd)
         izd = ixyz(3,ll+ipd)
*
      Do 200 ipc = 1, nElem(lc)
         ixc = ixyz(1,kk+ipc)
         iyc = ixyz(2,kk+ipc)
         izc = ixyz(3,kk+ipc)
*
         ixcd = ixc + ixd
         iycd = iyc + iyd
         izcd = izc + izd
*
      Do 300 ipb = 1, nElem(lb)
         ixb = ixyz(1,jj+ipb)
         iyb = ixyz(2,jj+ipb)
         izb = ixyz(3,jj+ipb)
*
         ixbcd = ixcd + ixb
         iybcd = iycd + iyb
         izbcd = izcd + izb
*
      Do 400 ipa = 1, nElem(la)
         ixa = ixyz(1,ii+ipa)
         iya = ixyz(2,ii+ipa)
         iza = ixyz(3,ii+ipa)
*
         ixabcd = ixbcd + ixa
         iyabcd = iybcd + iya
         izabcd = izbcd + iza
*
*
*        Compute all desired gradients with respect to an x-component.
*
         If (iyabcd.ne.0) Then
*
         If (nVec(1).eq.3) Then
           call ass3(xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &               xyz2D0(1,1,iza,izb,izc,izd,3),
     &               xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)),
     &               xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),
     &               xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(3,1)),
     &               PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),
     &               Temp(Ind2(2,1)),Temp(Ind2(3,1)),nT,nRys)
         Else If (nVec(1).eq.2) Then
           call ass2(xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &               xyz2D0(1,1,iza,izb,izc,izd,3),
     &               xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)),
     &               xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),
     &               PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),
     &               Temp(Ind2(2,1)),nT,nRys)
         Else If (nVec(1).eq.1) Then
           call ass1(xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &               xyz2D0(1,1,iza,izb,izc,izd,3),
     &               xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)),
     &               PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),
     &               nT,nRys)
         End If
*
         Else
*
         If (nVec(1).eq.3) Then
           call ass3a(xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)),
     &                xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),
     &                xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(3,1)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),
     &                Temp(Ind2(2,1)),Temp(Ind2(3,1)),nT,nRys)
         Else If (nVec(1).eq.2) Then
           call ass2a(xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)),
     &                xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),
     &                Temp(Ind2(2,1)),nT,nRys)
         Else If (nVec(1).eq.1) Then
           call ass1a(xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),
     &                nT,nRys)
         End If
*
         End If
*
*        Compute all desired gradients with respect to an y-component.
*
         If (ixabcd.ne.0) Then
*
         If (nVec(2).eq.3) Then
            call ass3(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(3,2)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),
     &                Temp(Ind2(2,2)),Temp(Ind2(3,2)),nT,nRys)
         Else If (nVec(2).eq.2) Then
            call ass2(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),
     &                Temp(Ind2(2,2)),nT,nRys)
         Else If (nVec(2).eq.1) Then
            call ass1(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),
     &                nT,nRys)
         End If
*
         Else
*
         If (nVec(2).eq.3) Then
           call ass3a(xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(3,2)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),
     &                Temp(Ind2(2,2)),Temp(Ind2(3,2)),nT,nRys)
         Else If (nVec(2).eq.2) Then
           call ass2a(xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),
     &                Temp(Ind2(2,2)),nT,nRys)
         Else If (nVec(2).eq.1) Then
           call ass1a(xyz2D0(1,1,iza,izb,izc,izd,3),
     &                xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),
     &                nT,nRys)
         End If
*
         End If
*
*        Compute all desired gradients with respect to an z-component.
*
         If (ixabcd*iyabcd.ne.0) Then
*
         If (nVec(3).eq.3) Then
            call ass3(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                Temp(Ind2(2,3)),Temp(Ind2(3,3)),nT,nRys)
         Else If (nVec(3).eq.2) Then
            call ass2(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                Temp(Ind2(2,3)),nT,nRys)
         Else If (nVec(3).eq.1) Then
            call ass1(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                nT,nRys)
         End If
*
         Else If (ixabcd.eq.0.and.iyabcd.ne.0) Then
*
         If (nVec(3).eq.3) Then
           call ass3a(xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                Temp(Ind2(2,3)),Temp(Ind2(3,3)),nT,nRys)
         Else If (nVec(3).eq.2) Then
           call ass2a(xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                Temp(Ind2(2,3)),nT,nRys)
         Else If (nVec(3).eq.1) Then
           call ass1a(xyz2D0(1,1,iya,iyb,iyc,iyd,2),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                nT,nRys)
         End If
*
         Else If (iyabcd.eq.0.and.ixabcd.ne.0) Then
*
         If (nVec(3).eq.3) Then
           call ass3a(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                Temp(Ind2(2,3)),Temp(Ind2(3,3)),nT,nRys)
         Else If (nVec(3).eq.2) Then
           call ass2a(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                Temp(Ind2(2,3)),nT,nRys)
         Else If (nVec(3).eq.1) Then
           call ass1a(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),
     &                xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                nT,nRys)
         End If
*
         Else
*
         If (nVec(3).eq.3) Then
            call ass3b(xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                 xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                 xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),
     &                 PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                 Temp(Ind2(2,3)),Temp(Ind2(3,3)),nT,nRys)
         Else If (nVec(3).eq.2) Then
            call ass2b(xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                 xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),
     &                 PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                 Temp(Ind2(2,3)),nT,nRys)
         Else If (nVec(3).eq.1) Then
            call ass1b(xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),
     &                 PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),
     &                 nT,nRys)
         End If
*
         End If
*
 400  Continue
*
 300  Continue
*
 200  Continue
*
 100  Continue
*
      Return
      End
