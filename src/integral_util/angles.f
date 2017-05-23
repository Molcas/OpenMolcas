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
      SubRoutine Angles(Lbls,xyz,mCentr,rtrnc,Max_Center)
************************************************************************
*                                                                      *
* Object: to compute angles from a list of coordinates.                *
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
#include "Molcas.fh"
      Real*8 xyz(3,mCentr)
      Character*(LENIN) Lbls(mCentr)
      Logical Type
*
      Lu=6
      iRout = 126
      iPrint = nPrint(iRout)
      If (mCentr.gt.Max_Center) Go To 99
*
      Type = .False.
*-----The center atom
      Do 52 ic = 1, mCentr
         x1 = xyz(1,ic)
         y1 = xyz(2,ic)
         z1 = xyz(3,ic)
         Do 53 jc = 1, mCentr
            If (jc.eq.ic) Go To 53
            x2 = xyz(1,jc)
            y2 = xyz(2,jc)
            z2 = xyz(3,jc)
            r1 = Sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            If (r1.gt.rtrnc.or.r1.eq.Zero) Go To 53
            Do 54 kc = jc+1, mCentr
               If (kc.eq.ic) Go To 54
               x3 = xyz(1,kc)
               y3 = xyz(2,kc)
               z3 = xyz(3,kc)
               r2 = Sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
               If (r2.gt.rtrnc.or.r2.eq.Zero) Go To 54
               arg=  ((x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)+
     &                   (z2-z1)*(z3-z1))/(r1*r2)
               If (abs(arg).gt.One) arg=sign(One,arg)
               Phi= 180.D0 * ACos(arg) / Pi
               If (.Not.Type) Then
                  Type = .True.
                  Write (Lu,*)
                  Write (Lu,'(19X,A)')
     &                  ' ************************************** '
                  Write (Lu,'(19X,A)')
     &                  ' *    Valence Bond Angles / Degree    * '
                  Write (Lu,'(19X,A)')
     &                  ' ************************************** '
                  Write (Lu,'(19X,A)')
     &                   '       Atom centers                 Phi'
               End If
               Write (Lu,'(21X,3(I2,1X,A,2X),1X,F6.2)')
     &               jc,Lbls(jc),ic,Lbls(ic),kc,Lbls(kc),Phi
 54         Continue
 53      Continue
 52   Continue
*
 99   Continue
      Return
      End
