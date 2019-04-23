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
      Subroutine Tranca(Cavxyz,Cavsph,lMax,CarSph)
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Dimension Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6),
     &          Cavsph( (lMax+1)**2 )
      Logical CarSph
*
      iOff1 = 1
      iOff2 = 1
      Do 100 m = 0, lMax
         nElem = (m+1)*(m+2)/2
*        Call RecPrt('Car-->Sph',' ',RSph(ipSph(m)),nElem,nElem)
         If (CarSph) Then
            call dcopy_(2*m+1,[Zero],0,Cavsph(iOff2),1)
*           Call RecPrt('Cartesian',' ',Cavxyz(iOff1),1,nElem)
            Call dGeMV_('T',nElem,2*m+1,
     &                 One,RSph(ipSph(m)),nElem,
     &                     Cavxyz(iOff1),1,
     &                 Zero,CavSph(iOff2),1)
*           Call RecPrt('Spherical',' ',Cavsph(iOff2),1,2*m+1)
         Else
            call dcopy_(nElem,[Zero],0,Cavxyz(iOff1),1)
*           Call RecPrt('Spherical',' ',Cavsph(iOff2),1,2*m+1)
            Call dGeMV_('N',nElem,2*m+1,
     &                 One,RSph(ipSph(m)),nElem,
     &                     Cavsph(iOff2),1,
     &                 Zero,Cavxyz(iOff1),1)
*           Call RecPrt('Cartesian',' ',Cavxyz(iOff1),1,nElem)
         End If
         iOff1 = iOff1 + nElem
         iOff2 = iOff2 + 2*m+1
 100  Continue
*
      Return
      End
