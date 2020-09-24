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
      SubRoutine HRR2Da(Arr1,nVec,nabMax,ncdMax,Arr2,A,B,la,lb,lc,ld,
     &                 IfGrad)
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
*             Modified to recurrence algorithm, February '92.          *
*             Improved algorithm, June '92.                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 A(3), B(3), Arr1(nVec,3,0:nabMax,0:ncdMax),
     &                   Arr2(nVec,0:la+1,0:lb+1,0:ncdMax,3)
      Logical IfGrad(3,4)
*
      iRout = 233
      iPrint = nPrint(iRout)
*     iQ = 0
*
      Do 10 iCar = 1, 3
         lla = 0
         If(IfGrad(iCar,1)) lla = 1
         llb = 0
         If(IfGrad(iCar,2)) llb = 1
         llcd = 0
         If(IfGrad(iCar,3).or.IfGrad(iCar,4)) llcd = 1
*
         AB = A(iCar)-B(iCar)
         If (AB.eq.Zero) Then
            Do 100 icd = 0, lc+ld+llcd
*-----------Using the identity
               Do 200 ia = 0, la+lla
                  Do 210 ib = 0, lb+llb
                  iab = ia + ib
                  If (iab.gt.la+lb+Max(lla,llb)) Go To 210
                  do i=1,nVec
                  Arr2(i,ia,ib,icd,iCar)=Arr1(i,iCar,iab,icd)
                  enddo
210               Continue
200            Continue
100         Continue
         Else
            If (la.ge.lb) Then
               Do 101 icd = 0, lc+ld+llcd
*-----------------Move the first row I(ia,0)
                  Do 20 ia = 0, la+lb+Max(lla,llb)
                     ja = ia
                     jb = 0
                     If (ja.gt.la+1) Then
                        ja = ja - (la+2)
                        jb = 1
                     End If
                     do i=1,nVec
                     Arr2(i,ja,jb,icd,iCar)=Arr1(i,iCar,ia,icd)
                     end do
 20               Continue
*-----------------Generate I(ia,ib) for fixed ib
                  Do 30 ib = 1, lb + llb
                     Do 31 ia = la+lb+Max(lla,llb)-ib, 0, -1
                        ja = ia
                        jb = ib
                        mb = ib-1
                        If (ja.gt.la+1) Then
                           ja = ja - (la+2)
                           jb = jb + 1
                           mb = mb + 1
                        End If
                        ma = ja
                        ka = ia+1
                        kb = ib-1
                        If (ka.gt.la+1) Then
                           ka = ka - (la+2)
                           kb = kb + 1
                        End If
                        Call DZaXpY(nVec,AB,Arr2(1,ma,mb,icd,iCar),1,
     &                                      Arr2(1,ka,kb,icd,iCar),1,
     &                                      Arr2(1,ja,jb,icd,iCar),1)
 31                  Continue
 30               Continue
 101           Continue
            Else
               AB = -AB
               Do 102 icd = 0, lc+ld+llcd
*-----------------Move the first row I(0,ib)
                  Do 40 ib = 0, la+lb+Max(lla,llb)
                     jb = ib
                     ja = 0
                     If (jb.gt.lb+1) Then
                        jb = jb - (lb+2)
                        ja = 1
                     End If
                     do i=1,nVec
                     Arr2(i,ja,jb,icd,iCar)=Arr1(i,iCar,ib,icd)
                     end do
 40               Continue
*-----------------Generate I(ia,ib) for fixed ia
                  Do 50 ia = 1, la + lla
                     Do 51 ib = la+lb+Max(lla,llb)-ia, 0, -1
                        jb = ib
                        ja = ia
                        ma = ia-1
                        If (jb.gt.lb+1) Then
                           jb = jb - (lb+2)
                           ja = ja + 1
                           ma = ma + 1
                        End If
                        mb = jb
                        kb = ib+1
                        ka = ia-1
                        If (kb.gt.lb+1) Then
                           kb = kb - (lb+2)
                           ka = ka + 1
                        End If
                        Call DZaXpY(nVec,AB,Arr2(1,ma,mb,icd,iCar),1,
     &                                      Arr2(1,ka,kb,icd,iCar),1,
     &                                      Arr2(1,ja,jb,icd,iCar),1)
 51                  Continue
 50               Continue
 102           Continue
            End If
         End If
 10   Continue
*
      Return
      End
