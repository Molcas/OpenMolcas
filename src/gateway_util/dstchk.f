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
* Copyright (C) 2009, Roland Lindh                                     *
************************************************************************
      SubRoutine DstChk(xyz,Lbls,mCentr)
************************************************************************
*                                                                      *
* Object: to check that the structure is realistic.                    *
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
#include "angstr.fh"
#include "Molcas.fh"
      Real*8 xyz(3,mCentr)
      Character*(LENIN) Lbls(mCentr)
*
      lu=6
      iRout = 125
      iPrint = nPrint(iRout)
*
      If (mCentr.lt.5) Go To 99
*
      iLarge=0
      Do icc=1,mCentr
        if(index('1234567890',Lbls(icc)(2:2)).eq.0)iLarge=1
      enddo
      if(iLarge.eq.1) goto 99
*
      RMax=0.0D0
      RMin=1.0D10
      iiBST=0
      Do icc = 1, mCentr
         x1 = xyz(1,icc)
         y1 = xyz(2,icc)
         z1 = xyz(3,icc)
         Do jcc = 1, icc-1
            x2 = xyz(1,jcc)
            y2 = xyz(2,jcc)
            z2 = xyz(3,jcc)
            R = Sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            RMin=Min(RMin,R)
            RMax=Max(RMax,r)
         End Do
      End Do
*
      If (RMax*Angstr.lt.0.7D0) Then
         Write (lu,*)
     &     'All bonds shorter than 0.7 angstrom, this is probably '//
     &     'wrong!'
         Write (lu,*)
     &     'The program will stop execution. To proceed, correct the '
         Write (lu,*)
     &     'input or use the "Expert" keyword to force execution'
         Call AbEnd()
      End If
      If (RMin*Angstr.gt.2.8D0) Then
         Write (lu,*)
     &     'All bonds longer than 2.8 angstrom, this is probably wrong!'
         Write (lu,*)
     &     'The program will stop execution. To proceed, correct the '
         Write (lu,*)
     &     'input or use the "Expert" keyword to force execution'
         Call AbEnd()
      End If
*
 99   Continue
      Return
      End
