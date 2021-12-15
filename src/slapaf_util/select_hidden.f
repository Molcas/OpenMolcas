************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Select_Hidden(mTtAtm,nHidden,Coord,HiddenCoord,
     &                         iHiddenAN,nKept,rHidden,iPL)
      Implicit Real*8 (a-h,o-z)
*
*  Select among the hidden atoms the ones to be kept
*
#include "real.fh"
      Dimension Coord(3,mTtAtm),HiddenCoord(3,nHidden)
      Dimension iHiddenAN(nHidden)
*
*
*  Criteria: dMin < distance < rHidden
*
      dMax = rHidden ! bohr
      Do iHid = 1, nHidden
         X = HiddenCoord(1,iHid)
         Y = HiddenCoord(2,iHid)
         Z = HiddenCoord(3,iHid)
         iAN = iHiddenAN(iHid)
         iAtom = 0
10       iAtom = iAtom + 1
         Dist = sqrt((X-Coord(1,iAtom))**2+(Y-Coord(2,iAtom))**2+
     &               (Z-Coord(3,iAtom))**2)
         If(Dist .le. dMax) Then
            iHiddenAN(iHid) = - iAN
            nKept = nKept + 1
         End If
         If(iAtom.lt.mTtAtm .and. iHiddenAN(iHid).le.0) Goto 10
      End Do
*
*  The end
*
      If(iPL .gt. 3 .and. nKept .gt. 0) Write(6,'(A,i3,A)')
     &                ' Select_Hidden: ',nKept,' hidden atoms are kept'
*
      Return
      End
