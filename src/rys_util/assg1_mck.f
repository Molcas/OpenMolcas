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
      SubRoutine Assg1_mck(g1,nT,nRys,la,lb,lc,ld,xyz2D0,xyz2D1,
     &                 IfGrad,Index,mVec,Index2)
************************************************************************
*                                                                      *
* Object: to assemble the gradients of the ERI's.                      *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "iavec.fh"
      Real*8 g1(nT,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,
     &             (lc+1)*(lc+2)/2,(ld+1)*(ld+2)/2,9),
     &       xyz2D0(nRys,nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3),
     &       xyz2D1(nRys,nT,0:la  ,0:lb  ,0:lc  ,0:ld  ,9)
      Logical IfGrad(3,4)
      Integer Ind1(3), Ind2(3), Index(3,4),Index2(3,4)
*
*     Statement functions
*
      nElem(i) = (i+1)*(i+2)/2
*
      ka=(la+1)*(la+2)/2
      kb=(lb+1)*(lb+2)/2
      kc=(lc+1)*(lc+2)/2
      kd=(ld+1)*(ld+2)/2
      nG1=nT*9*ka*kb*kc*kd
      call dcopy_(nG1,[Zero],0,G1,1)
      Call ICOPY(12,[0],0,Index2,1)
*
      ii = la*(la+1)*(la+2)/6
      jj = lb*(lb+1)*(lb+2)/6
      kk = lc*(lc+1)*(lc+2)/6
      ll = ld*(ld+1)*(ld+2)/6
*
      Do 100 ipa = 1, nElem(la)
         ipaii=ipa+ii
         ixa = ixyz(1,ipaii)
         iya = ixyz(2,ipaii)
         iza = ixyz(3,ipaii)
*
      Do 200 ipb = 1, nElem(lb)
         ipbjj=ipb+jj
         ixb = ixyz(1,ipbjj)
         iyb = ixyz(2,ipbjj)
         izb = ixyz(3,ipbjj)
*
         ixab = ixa + ixb
         iyab = iya + iyb
         izab = iza + izb
*
      Do 300 ipc = 1, nElem(lc)
         ipckk=ipc+kk
         ixc = ixyz(1,ipckk)
         iyc = ixyz(2,ipckk)
         izc = ixyz(3,ipckk)
*
         ixabc = ixab + ixc
         iyabc = iyab + iyc
         izabc = izab + izc
*
      Do 400 ipd = 1, nElem(ld)
         ipdll=ipd+ll
         ixd = ixyz(1,ipdll)
         iyd = ixyz(2,ipdll)
         izd = ixyz(3,ipdll)
*
         ixabcd = ixabc + ixd
         iyabcd = iyabc + iyd
         izabcd = izabc + izd
*
*        Compute all desired gradients with respect to an x-component.
*
         mVec = 0
         nVec = 0
         Do 1000 iCent = 1, 4
            If (IfGrad(1,iCent)) Then
               mVec = mVec + 1
               nVec = nVec + 1
               Ind1(nVec) = 3*(Index(1,iCent)-1) + 1
               Ind2(nVec) = mVec
               Index2(1,iCent)=mVec
            End If
 1000    Continue
*
         If (iyabcd.ne.0) Then
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2) *
     &                  xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
                  tmp3 = tmp3 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1 = 0.0D0
               tmp2 = 0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2) *
     &                  xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1 = 0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2) *
     &                  xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
         End If
*
         Else
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
                  tmp3 = tmp3 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
         End If
*
         End If
*
*        Compute all desired gradients with respect to an y-component.
*
         nVec = 0
         Do 2000 iCent = 1, 4
            If (IfGrad(2,iCent)) Then
               mVec = mVec + 1
               nVec = nVec + 1
               Ind1(nVec) = 3*(Index(2,iCent)-1) + 2
               Ind2(nVec) = mVec
               Index2(2,iCent)=mVec
            End If
 2000    Continue
*
         If (ixabcd.ne.0) Then
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1) *
     &                  xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
                  tmp3 = tmp3 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1) *
     &                  xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1 = 0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1) *
     &                  xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
         End If
*
         Else
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
                  tmp3 = tmp3 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1 = 0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
         End If
*
         End If
*
*        Compute all desired gradients with respect to an z-component.
*
         nVec = 0
         Do 3000 iCent = 1, 4
            If (IfGrad(3,iCent)) Then
               mVec = mVec + 1
               nVec = nVec + 1
               Ind1(nVec) = 3*(Index(3,iCent)-1) + 3
               Ind2(nVec) =  mVec
               Index2(3,iCent)=mVec
            End If
 3000    Continue
*
         If (ixabcd*iyabcd.ne.0) Then
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1) *
     &                  xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1) *
     &                  xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1) *
     &                  xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
         End If
*
         Else If (ixabcd.eq.0.and.iyabcd.ne.0) Then
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
         End If
*
         Else If (iyabcd.eq.0.and.ixabcd.ne.0) Then
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1=0.0D0
               Do iRys = 1, nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)
                  tmp1 = tmp1 + tmp *
     &                  xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
         End If
*
         Else
*
         If (nVec.eq.3) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               tmp3=0.0D0
               Do iRys = 1, nRys
                  tmp1 = tmp1 + xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3 + xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
               g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) + tmp3
            End Do
         Else If (nVec.eq.2) Then
            Do iT = 1, nT
               tmp1=0.0D0
               tmp2=0.0D0
               Do iRys = 1, nRys
                  tmp1 = tmp1 + xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2 + xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
               g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) + tmp2
            End Do
         Else If (nVec.eq.1) Then
            Do iT = 1, nT
               tmp1=0.0D0
               Do iRys = 1, nRys
                  tmp1 = tmp1 + xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
               End Do
               g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) =
     &                    g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) + tmp1
            End Do
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
