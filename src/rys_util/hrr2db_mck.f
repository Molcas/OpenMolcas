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
* Copyright (C) 1991,1992, Roland Lindh                                *
************************************************************************
      SubRoutine Hrr2Db_mck(Arr1,nVec,ncdMax,
     &                 Arr2,C,D,la,lb,lc,ld,IfHss,
     &                 IfGrd ,nt,nrys)
************************************************************************
*                                                                      *
* Object: to apply the transfer equation to the 2D-integrals.          *
*         The transformation is in place and the recursion             *
*         is replaced with the indentity when applicable.              *
*                                                                      *
* Called from: HrrCtl                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy   (ESSL)                                          *
*              DZaXpY  (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September '91                                            *
*             Modified to recurrence algorithm, February '92           *
*             Improved algorithm, June '92.                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
#include "real.fh"
      Real*8 C(3), D(3),
     &                   Arr1(nVec,0:la+2,0:lb+2,0:ncdMax,3),
     &                   Arr2(nVec,0:la+2,0:lb+2,0:lc+2,0:ld+2,3)
      Logical IfHss(4,3,4,3),ifgrd(3,4)
*
c     iRout = 233
c     iPrint = nPrint(iRout)
*
      Do 10 iCar = 1, 3
         llc = 0
         lld = 0
         lla = 0
         llb = 0
         If(IfGrd(iCar,3)) llc = Max(llc,1)
         If(IfHss(3,iCar,3,iCar)) llc = 2
         If(IfGrd(iCar,4)) lld = Max(1,lld)
         If(IfHss(4,iCar,4,iCar)) lld = 2
         If(IfGrd(iCar,1)) lla = Max(1,lla)
         If(IfHss(1,iCar,1,iCar)) lla = 2
         If(IfGrd(iCar,2)) llb = Max(llb,1)
         If(IfHss(2,iCar,2,iCar)) llb = 2
         CD = C(iCar)-D(iCar)
         If (CD.eq.Zero) Then
            Do 100 ia = 0, la+lla
            Do 100 ib = 0, lb+llb
               If (ia+ib.gt.la+lb+Max(lla,llb)) Go To 100
*--------------Using the identity
               Do 200 ic = 0, lc+llc
                  Do 210 id = 0, ld+lld
                  icd = ic + id
                  If (icd.gt.lc+ld+Max(llc,lld)) Go To 210
                  do i=1,nVec
                  Arr2(i,ia,ib,ic,id,iCar)=Arr1(i,ia,ib,icd,iCar)
                  enddo
210               Continue
200            Continue
100         Continue
         Else
            If (lc.ge.ld) Then
               Do 101 ia = 0, la+lla
               Do 101 ib = 0, lb+llb
                  If (ia+ib.gt.la+lb+Max(lla,llb)) Go To 101
*-----------------Move the first row I(ic,0)
                  Do 20 ic = 0, lc+ld+Max(llc,lld)
                     jc = ic
                     jd = 0
                     If (jc.gt.lc+2) Then
                        jc = jc - (lc+3)
                        jd = 1
                     End If
                     do i=1,nVec
                     Arr2(i,ia,ib,jc,jd,iCar)=Arr1(i,ia,ib,ic,iCar)
                     enddo
 20               Continue
*-----------------Generate I(ic,id) for fixed id
                  Do 30 id = 1, ld + lld
                     Do 31 ic = lc+ld+Max(llc,lld)-id, 0, -1
                        jc = ic
                        jd = id
                        md = id-1
                        If (jc.gt.lc+2) Then
                           jc = jc - (lc+3)
                           jd = jd + 1
                           md = md + 1
                        End If
                        mc = jc
                        kc = ic+1
                        kd = id-1
                        If (kc.gt.lc+2) Then
                           kc = kc - (lc+3)
                           kd = kd + 1
                        End If
                        Call DZaXpY(nVec,CD,Arr2(1,ia,ib,mc,md,iCar),1,
     &                                      Arr2(1,ia,ib,kc,kd,iCar),1,
     &                                      Arr2(1,ia,ib,jc,jd,iCar),1)
 31                  Continue
 30               Continue
 101           Continue
            Else
               CD = -CD
               Do 102 ia = 0, la+lla
               Do 102 ib = 0, lb+llb
                  If (ia+ib.gt.la+lb+Max(lla,llb)) Go To 102
*-----------------Move the first row I(0,id)
                  Do 40 id = 0, lc+ld+Max(llc,lld)
                     jd = id
                     jc = 0
                     If (jd.gt.ld+2) Then
                        jd = jd - (ld+3)
                        jc = 1
                     End If
                     do i=1,nVec
                     Arr2(i,ia,ib,jc,jd,iCar)=Arr1(i,ia,ib,id,iCar)
                     enddo
 40               Continue
*-----------------Generate I(ic,id) for fixed ic
                  Do 50 ic = 1, lc + llc
                     Do 51 id = lc+ld+Max(llc,lld)-ic, 0, -1
                        jd = id
                        jc = ic
                        mc = ic-1
                        If (jd.gt.ld+2) Then
                           jd = jd - (ld+3)
                           jc = jc + 1
                           mc = mc + 1
                        End If
                        md = jd
                        kd = id+1
                        kc = ic-1
                        If (kd.gt.ld+2) Then
                           kd = kd - (ld+3)
                           kc = kc + 1
                        End If
                        Call DZaXpY(nVec,CD,Arr2(1,ia,ib,mc,md,iCar),1,
     &                                      Arr2(1,ia,ib,kc,kd,iCar),1,
     &                                      Arr2(1,ia,ib,jc,jd,iCar),1)
 51                  Continue
 50               Continue
 102           Continue
            End If
         End If
 10   Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nt)
         Call Unused_integer(nrys)
      End If
      End
