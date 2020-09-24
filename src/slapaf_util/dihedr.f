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
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      SubRoutine Dihedr(Lbls,xyz,mCentr,rtrnc,Max_Center)
************************************************************************
*                                                                      *
* Object: to compute dihedral angles from a list of coordinates.       *
*                                                                      *
* Called from: Input                                                   *
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
#include "Molcas.fh"
      Real*8 xyz(3,mCentr), Bt(3,4), Coor(3,4)
      Character*(LENIN) Lbls(mCentr)
      Character*8 Label
      Logical Type
      Dimension Dummy(1)
*
      Lu=6
      iRout = 127
      iPrint = nPrint(iRout)
      Label=' '
      If (mCentr.gt.Max_Center) Go To 99
*
      Thr = 1.0D-12
      Type = .False.
      Do 452 ic = 1, mCentr
         x2 = xyz(1,ic)
         y2 = xyz(2,ic)
         z2 = xyz(3,ic)
         call dcopy_(3,xyz(1,ic),1,Coor(1,2),1)
         Do 453 jc = 1, mCentr
            If (jc.eq.ic) Go To 453
            x3 = xyz(1,jc)
            y3 = xyz(2,jc)
            z3 = xyz(3,jc)
            r2 = Sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
            If (r2.gt.rtrnc.or.r2.eq.Zero) Go To 453
*           Write (Lu,*)
            call dcopy_(3,xyz(1,jc),1,Coor(1,3),1)
            Do 454 kc = 1, mCentr
               If (kc.eq.ic) Go To 454
               If (kc.eq.jc) Go To 454
               x1 = xyz(1,kc)
               y1 = xyz(2,kc)
               z1 = xyz(3,kc)
               r1 = Sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
               If (r1.gt.rtrnc.or.r1.eq.Zero) Go To 454
               arg=  ((x1-x2)*(x3-x2)+(y1-y2)*(y3-y2)+
     &                   (z1-z2)*(z3-z2))/(r1*r2)
               If (abs(arg).gt.One) arg=sign(One,arg)
               If (One-abs(arg).lt.Thr) Go To 454
               Phi1= 180.D0 * ACos(arg) / Pi
               x12 = (y2-y1)*(z3-z2)-(y3-y2)*(z2-z1)
               y12 = (z2-z1)*(x3-x2)-(z3-z2)*(x2-x1)
               z12 = (x2-x1)*(y3-y2)-(x3-x2)*(y2-y1)
               r12 = Sqrt(x12**2+y12**2+z12**2)
               If (r12.eq.Zero) Go To 454
               call dcopy_(3,xyz(1,kc),1,Coor(1,1),1)
               Do 455 lc = kc+1, mCentr
                  If (lc.eq.ic) Go To 455
                  If (lc.eq.jc) Go To 455
                  If (lc.eq.kc) Go To 455
                  x4 = xyz(1,lc)
                  y4 = xyz(2,lc)
                  z4 = xyz(3,lc)
                  r3 = Sqrt((x4-x3)**2+(y4-y3)**2+(z4-z3)**2)
                  If (r3.gt.rtrnc.or.r3.eq.Zero) Go To 455
                  arg=  ((x2-x3)*(x4-x3)+(y2-y3)*(y4-y3)+
     &                      (z2-z3)*(z4-z3))/(r2*r3)
                  If (abs(arg).gt.One) arg=sign(One,arg)
                  If (One-abs(arg).lt.Thr) Go To 455
                  Phi2= 180.D0 * ACos(arg) / Pi
                  x23 = (y3-y2)*(z4-z3)-(y4-y3)*(z3-z2)
                  y23 = (z3-z2)*(x4-x3)-(z4-z3)*(x3-x2)
                  z23 = (x3-x2)*(y4-y3)-(x4-x3)*(y3-y2)
                  r23 = Sqrt(x23**2+y23**2+z23**2)
                  If (r23.eq.Zero) Go To 455
                  call dcopy_(3,xyz(1,lc),1,Coor(1,4),1)
C                 arg=  (x12*x23+y12*y23+z12*z23)/
C    &                                     (r12*r23)
C                 If (abs(arg).gt.One) arg=sign(One,arg)
C                 Phi12= 180.D0 * ACos(arg) / Pi
                  Call Trsn(Coor,4,Tau,Bt,.False.,.False.,Label,
     &                      Dummy,.False.)
                  Phi12 = 180.0D+00*Tau/Pi
                  If (.Not.Type) Then
                     Type = .True.
                     Write (Lu,*)
                     Write (Lu,'(10X,A)')
     &' ***************************************************************'
                     Write (Lu,'(10X,A)')
     &' *              Valence Dihedral Angles / Degree               *'
                     Write (Lu,'(10X,A)')
     &' ***************************************************************'
                     Write (Lu,'(7X,A,A)')
     &                     '             Atom centers       ',
     &                     '                Phi1     Phi2     Theta '
                  End If
                  Write (Lu,'(10X,4(I2,1X,A,2X),1X,3(F7.2,2X))')
     &                  kc,Lbls(kc),
     &                  ic,Lbls(ic),
     &                  jc,Lbls(jc),
     &                  lc,Lbls(lc),
     &               Phi1, Phi2, Phi12
 455           Continue
 454         Continue
 453      Continue
 452   Continue
*
 99   Continue
      Return
      End
