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
      Subroutine PPMmG(nHer,MmPPG,la,lb,lr)
*
      nElem(i) = (i+1)*(i+2)/2
*
      nHer=0
      MmPPG=0
*
      laplb=Max(nElem(la+1),nElem(lb))**2
      MmPPG=MmPPG+2*laplb
*
      If (la.gt.0) Then
         lamlb=Max(nElem(la-1),nElem(lb))**2
      Else
         lamlb=0
      End If
      MmPPG=MmPPG+2*lamlb
*
      lalbp=Max(nElem(la),nElem(lb+1))**2
      MmPPG=MmPPG+2*laplb
*
      If (lb.gt.0) Then
         lalbm=Max(nElem(la),nElem(lb-1))**2
      Else
         lalbm=0
      End If
      MmPPG=MmPPG+2*lalbm
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
