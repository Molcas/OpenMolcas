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
      SubRoutine Rys2Dg(xyz2D0,nT,nRys,la,lb,lc,ld,xyz2D1,IfGrad,
     &                  IndGrd,Coora,Alpha,Beta,Gamma,Delta,nZeta,
     &                  nEta,Scrtch,Temp,Index,ExpX,ExpY,mZeta,mEta)
************************************************************************
*                                                                      *
* Object: to compute the gradients of the 2D-integrals.                *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External ExpX, ExpY
#include "print.fh"
#include "real.fh"
      Real*8 xyz2D0(nRys*nT,0:la+1,0:lb+1,0:lc+1,0:ld+1,3),
     &       xyz2D1(nRys*nT,0:la  ,0:lb  ,0:lc  ,0:ld  ,3,3),
     &       Coora(3,4),
     &       Alpha(nZeta), Beta(nZeta), Gamma(nEta), Delta(nEta),
     &       Scrtch(nRys*nT), Temp(nT)
      Logical IfGrad(3,4), EQ
      Integer IndGrd(3,4), Ind1(3), Ind2(3), Index(3,4)
*
#ifdef _DEBUGPRINT_
      iRout = 249
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In Rys2Dg: Alpha',' ',Alpha,1,nZeta)
         Call RecPrt(' In Rys2Dg: Beta ',' ',Beta ,1,nZeta)
         Call RecPrt(' In Rys2Dg: Gamma',' ',Gamma,1,nEta )
         Call RecPrt(' In Rys2Dg: Delta',' ',Delta,1,nEta )
         Write (6,*) ' IfGrad=',IfGrad
         Write (6,*) ' IndGrd=',IndGrd
      End If
#endif
      tOne = -One
      tTwo = Two
      nx = 0
      ny = 0
      nz = 0
      Call ICopy(12,[0],0,Index,1)
*
*     Differentiate with respect to the first center
*
      If (IfGrad(1,1) .or.IfGrad(2,1) .or.
     &                    IfGrad(3,1)) Then
         Call ExpX(Temp  ,mZeta,mEta,Alpha,One)
         Call Exp_2(Scrtch,nRys,nT,Temp,One)
*        If (iPrint.ge.99) Call RecPrt(
*    &      'Expanded exponents (alpha)',' ',Scrtch,
*    &      nT,nRys)
      End If
      nVec = 0
      If (IfGrad(1,1)) Then
         nx = nx + 1
         nVec = nVec + 1
         Ind1(nVec) = nx
         Ind2(nVec) = 1
         Index(1,1) = nx
      End If
      If (IfGrad(2,1)) Then
         ny = ny + 1
         nVec = nVec + 1
         Ind1(nVec) = ny
         Ind2(nVec) = 2
         Index(2,1) = ny
      End If
      If (IfGrad(3,1)) Then
         nz = nz + 1
         nVec = nVec + 1
         Ind1(nVec) = nz
         Ind2(nVec) = 3
         Index(3,1) = nz
      End If
      If (nVec.eq.0) Go To 211
*
      Do 101 id = 0, ld
         Do 201 ic = 0, lc
            Do 301 ib = 0, lb
               If (nVec.eq.3) Then
                  Do 501 iVec = 1, nT*nRys
                     xyz2D1(iVec,0,ib,ic,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,1,ib,ic,id,Ind2(1))
                     xyz2D1(iVec,0,ib,ic,id,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,1,ib,ic,id,Ind2(2))
                     xyz2D1(iVec,0,ib,ic,id,Ind2(3),Ind1(3)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,1,ib,ic,id,Ind2(3))
 501              Continue
                  If (la.ge.1) Then
                     Fact = tOne
                     Do 511 ia = 1, la
                        Do 521 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia+1,ib,ic,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia-1,ib,ic,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia+1,ib,ic,id,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia-1,ib,ic,id,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia+1,ib,ic,id,Ind2(3))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia-1,ib,ic,id,Ind2(3))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(3),Ind1(3)) =
     &                          tmp1 + tmp2
 521                    Continue
                        Fact = Fact + tOne
 511                 Continue
                  End If
               Else If (nVec.eq.2) Then
                  Do 601 iVec = 1, nT*nRys
                     xyz2D1(iVec,0,ib,ic,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,1,ib,ic,id,Ind2(1))
                     xyz2D1(iVec,0,ib,ic,id,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,1,ib,ic,id,Ind2(2))
 601              Continue
                  If (la.ge.1) Then
                     Fact = tOne
                     Do 611 ia = 1, la
                        Do 621 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia+1,ib,ic,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia-1,ib,ic,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia+1,ib,ic,id,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia-1,ib,ic,id,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
 621                    Continue
                        Fact = Fact + tOne
 611                 Continue
                  End If
               Else If (nVec.eq.1) Then
                  Do 701 iVec = 1, nT*nRys
                     xyz2D1(iVec,0,ib,ic,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,1,ib,ic,id,Ind2(1))
 701              Continue
                  If (la.ge.1) Then
                     Fact = tOne
                     Do 711 ia = 1, la
                        Do 721 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia+1,ib,ic,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia-1,ib,ic,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
 721                    Continue
                        Fact = Fact + tOne
 711                 Continue
                  End If
               End If
 301        Continue
 201     Continue
 101  Continue
*
*     Differentiate with respect to the second center
*
 211  Continue
      If (IfGrad(1,2) .or.IfGrad(2,2) .or.
     &                    IfGrad(3,2)) Then
         Call ExpX(Temp  ,mZeta,mEta,Beta,One)
         Call Exp_2(Scrtch,nRys,nT,Temp,One)
*        If (iPrint.ge.99) Call RecPrt(
*    &      'Expanded exponents (beta) ',' ',Scrtch,
*    &      nT,nRys)
      End If
      nVec = 0
      If (IfGrad(1,2)) Then
         nx = nx + 1
         nVec = nVec + 1
         Ind1(nVec) = nx
         Ind2(nVec) = 1
         Index(1,2) = nx
      End If
      If (IfGrad(2,2)) Then
         ny = ny + 1
         nVec = nVec + 1
         Ind1(nVec) = ny
         Ind2(nVec) = 2
         Index(2,2) = ny
      End If
      If (IfGrad(3,2)) Then
         nz = nz + 1
         nVec = nVec + 1
         Ind1(nVec) = nz
         Ind2(nVec) = 3
         Index(3,2) = nz
      End If
      If (nVec.eq.0) Go To 311
*
      Do 102 id = 0, ld
         Do 202 ic = 0, lc
            Do 302 ia = 0, la
               If (nVec.eq.3) Then
                  Do 502 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,0,ic,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,1,ic,id,Ind2(1))
                     xyz2D1(iVec,ia,0,ic,id,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,1,ic,id,Ind2(2))
                     xyz2D1(iVec,ia,0,ic,id,Ind2(3),Ind1(3)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,1,ic,id,Ind2(3))
 502              Continue
                  If (lb.ge.1) Then
                     Fact = tOne
                     Do 512 ib = 1, lb
                        Do 522 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib+1,ic,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib-1,ic,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib+1,ic,id,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib-1,ic,id,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib+1,ic,id,Ind2(3))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib-1,ic,id,Ind2(3))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(3),Ind1(3)) =
     &                          tmp1 + tmp2
 522                    Continue
                        Fact = Fact + tOne
 512                 Continue
                  End If
               Else If (nVec.eq.2) Then
                  Do 602 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,0,ic,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,1,ic,id,Ind2(1))
                     xyz2D1(iVec,ia,0,ic,id,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,1,ic,id,Ind2(2))
 602              Continue
                  If (lb.ge.1) Then
                     Fact = tOne
                     Do 612 ib = 1, lb
                        Do 622 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib+1,ic,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib-1,ic,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib+1,ic,id,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib-1,ic,id,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
 622                    Continue
                        Fact = Fact + tOne
 612                 Continue
                  End If
               Else If (nVec.eq.1) Then
                  Do 702 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,0,ic,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,1,ic,id,Ind2(1))
 702              Continue
                  If (lb.ge.1) Then
                     Fact = tOne
                     Do 712 ib = 1, lb
                        Do 722 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib+1,ic,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib-1,ic,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
 722                    Continue
                        Fact = Fact + tOne
 712                 Continue
                  End If
               End If
 302        Continue
 202     Continue
 102  Continue
*
*     Differentiate with respect to the third center
*
 311  Continue
      If (IfGrad(1,3) .or.IfGrad(2,3) .or.
     &                    IfGrad(3,3)) Then
         Call ExpY(Temp  ,mZeta,mEta,Gamma,One)
         Call Exp_2(Scrtch,nRys,nT,Temp,One)
*        If (iPrint.ge.99) Call RecPrt(
*    &      'Expanded exponents (gamma)',' ',Scrtch,
*    &      nT,nRys)
      End If
      nVec = 0
      If (IfGrad(1,3)) Then
         nx = nx + 1
         nVec = nVec + 1
         Ind1(nVec) = nx
         Ind2(nVec) = 1
         Index(1,3) = nx
      End If
      If (IfGrad(2,3)) Then
         ny = ny + 1
         nVec = nVec + 1
         Ind1(nVec) = ny
         Ind2(nVec) = 2
         Index(2,3) = ny
      End If
      If (IfGrad(3,3)) Then
         nz = nz + 1
         nVec = nVec + 1
         Ind1(nVec) = nz
         Ind2(nVec) = 3
         Index(3,3) = nz
      End If
      If (nVec.eq.0) Go To 411
*
      Do 103 id = 0, ld
         Do 203 ib = 0, lb
            Do 303 ia = 0, la
               If (nVec.eq.3) Then
                  Do 503 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,ib,0,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,1,id,Ind2(1))
                     xyz2D1(iVec,ia,ib,0,id,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,1,id,Ind2(2))
                     xyz2D1(iVec,ia,ib,0,id,Ind2(3),Ind1(3)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,1,id,Ind2(3))
 503              Continue
                  If (lc.ge.1) Then
                     Fact = tOne
                     Do 513 ic = 1, lc
                        Do 523 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic+1,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic-1,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic+1,id,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic-1,id,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic+1,id,Ind2(3))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic-1,id,Ind2(3))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(3),Ind1(3)) =
     &                          tmp1 + tmp2
 523                    Continue
                        Fact = Fact + tOne
 513                 Continue
                  End If
               Else If (nVec.eq.2) Then
                  Do 603 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,ib,0,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,1,id,Ind2(1))
                     xyz2D1(iVec,ia,ib,0,id,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,1,id,Ind2(2))
 603              Continue
                  If (lc.ge.1) Then
                     Fact = tOne
                     Do 613 ic = 1, lc
                        Do 623 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic+1,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic-1,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic+1,id,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic-1,id,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
 623                    Continue
                        Fact = Fact + tOne
 613                 Continue
                  End If
               Else If (nVec.eq.1) Then
                  Do 703 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,ib,0,id,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,1,id,Ind2(1))
 703              Continue
                  If (lc.ge.1) Then
                     Fact = tOne
                     Do 713 ic = 1, lc
                        Do 723 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic+1,id,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic-1,id,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
 723                    Continue
                        Fact = Fact + tOne
 713                 Continue
                  End If
               End If
 303        Continue
 203     Continue
 103  Continue
*
*     Differentiate with respect to the fourth center
*
 411  Continue
      If (IfGrad(1,4) .or. IfGrad(2,4) .or.
     &                     IfGrad(3,4)) Then
         Call ExpY(Temp  ,mZeta,mEta,Delta,One)
         Call Exp_2(Scrtch,nRys,nT,Temp,One)
*        If (iPrint.ge.99) Call RecPrt(
*    &      'Expanded exponents (delta)',' ',Scrtch,
*    &      nT,nRys)
      End If
      nVec = 0
      If (IfGrad(1,4)) Then
         nx = nx + 1
         nVec = nVec + 1
         Ind1(nVec) = nx
         Ind2(nVec) = 1
         Index(1,4) = nx
      End If
      If (IfGrad(2,4)) Then
         ny = ny + 1
         nVec = nVec + 1
         Ind1(nVec) = ny
         Ind2(nVec) = 2
         Index(2,4) = ny
      End If
      If (IfGrad(3,4)) Then
         nz = nz + 1
         nVec = nVec + 1
         Ind1(nVec) = nz
         Ind2(nVec) = 3
         Index(3,4) = nz
      End If
      If (nVec.eq.0) Go To 999
*
      Do 104 ic = 0, lc
         Do 204 ib = 0, lb
            Do 304 ia = 0, la
               If (nVec.eq.3) Then
                  Do 504 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,ib,ic,0,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,1,Ind2(1))
                     xyz2D1(iVec,ia,ib,ic,0,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,1,Ind2(2))
                     xyz2D1(iVec,ia,ib,ic,0,Ind2(3),Ind1(3)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,1,Ind2(3))
 504              Continue
                  If (ld.ge.1) Then
                     Fact = tOne
                     Do 514 id = 1, ld
                        Do 524 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic,id+1,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic,id-1,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic,id+1,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic,id-1,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic,id+1,Ind2(3))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic,id-1,Ind2(3))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(3),Ind1(3)) =
     &                          tmp1 + tmp2
 524                    Continue
                        Fact = Fact + tOne
 514                 Continue
                  End If
               Else If (nVec.eq.2) Then
                  Do 604 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,ib,ic,0,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,1,Ind2(1))
                     xyz2D1(iVec,ia,ib,ic,0,Ind2(2),Ind1(2)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,1,Ind2(2))
 604              Continue
                  If (ld.ge.1) Then
                     Fact = tOne
                     Do 614 id = 1, ld
                        Do 624 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic,id+1,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic,id-1,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic,id+1,Ind2(2))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic,id-1,Ind2(2))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(2),Ind1(2)) =
     &                          tmp1 + tmp2
 624                    Continue
                        Fact = Fact + tOne
 614                 Continue
                  End If
               Else If (nVec.eq.1) Then
                  Do 704 iVec = 1, nT*nRys
                     xyz2D1(iVec,ia,ib,ic,0,Ind2(1),Ind1(1)) =
     &                   tTwo * Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,1,Ind2(1))
 704              Continue
                  If (ld.ge.1) Then
                     Fact = tOne
                     Do 714 id = 1, ld
                        Do 724 iVec = 1, nT*nRys
                           tmp1 = tTwo * Scrtch(iVec) *
     &                          xyz2D0(iVec,ia,ib,ic,id+1,Ind2(1))
                           tmp2 = Fact *
     &                          xyz2D0(iVec,ia,ib,ic,id-1,Ind2(1))
                           xyz2D1(iVec,ia,ib,ic,id,Ind2(1),Ind1(1)) =
     &                          tmp1 + tmp2
 724                    Continue
                        Fact = Fact + tOne
 714                 Continue
                  End If
               End If
 304        Continue
 204     Continue
 104  Continue
*
 999  Continue
*
*-----Sum over common centers
*
      Do 1000 iCent = 1, 3
         Do 2000 jCent = iCent+1, 4
            If (EQ(Coora(1,iCent),Coora(1,jCent))) Then
               Do 3000 iCar = 1, 3
*
                  If (IfGrad(iCar,iCent) .and.
     &                IfGrad(iCar,jCent) ) Then
*--------------------Change flags so gradient will not be assembled and
*                    that there will be no contribution to the gradient.
                     IfGrad(iCar,jCent) = .False.
                     IndGrd(iCar,jCent) = 0
                     i1 = Index(iCar,iCent)
                     i2 = Index(iCar,jCent)
                     Call DaXpY_(nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),
     &                          One,
     &                          xyz2D1(1,0,0,0,0,iCar,i2),1,
     &                          xyz2D1(1,0,0,0,0,iCar,i1),1)
                  End If
*
 3000          Continue
            End If
 2000    Continue
 1000 Continue
*
c     If (iPrint.ge.49) Then
c        Do 900 iCn = 1, 4
c           Do 901 iCar = 1, 3
c              If (IfGrad(iCar,iCn)) Then
c                 ij = Index(iCar,iCn)
c                 Do 902 ia = 0, la
c                    Do 903 ib = 0, lb
c                       Do 904 ic = 0, lc
c                          Do 905 id = 0, ld
c                             Write
c    &                         (Label,'(A,4(I2,'',''),A,'','',I2,A)')
c    &                           ' xyz2D1(',
c    &                           ia,ib,ic,id,ch(iCar),iCn,')'
c                             If (iPrint.ge.99) Then
c                                Call RecPrt(Label,' ',xyz2d1(1,ia,ib,
c    &                                    ic,id,iCar,ij),nT,nRys)
c                             Else
c                                Write (*,'(A)') Label
c                                Write (*,*) DDot_(nT*nRys,
c    &                            xyz2d1(1,ia,ib,ic,id,iCar,ij),1,
c    &                            xyz2d1(1,ia,ib,ic,id,iCar,ij),1)
c                             End If
c905                       Continue
c904                    Continue
c903                 Continue
c902              Continue
c              End If
c901        Continue
c900     Continue
c     End If
      Return
      End
