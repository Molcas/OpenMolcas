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
      Subroutine Compute_M(ZA,nAtoms,RA,Z_Tot,T,M)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 ZA(nAtoms), RA(3,nAtoms), T(3), M(3,3)
*                                                                      *
************************************************************************
*                                                                      *
*---- Form the nuclear charge moment tensor
*
      Call FZero(M,9)
      Do iAtom = 1, nAtoms
         RTx=RA(1,iAtom)-T(1)
         RTy=RA(2,iAtom)-T(2)
         RTz=RA(3,iAtom)-T(3)
         M(1,1) = M(1,1) + ZA(iAtom) * (RTy**2+RTz**2)
         M(2,2) = M(2,2) + ZA(iAtom) * (RTx**2+RTz**2)
         M(3,3) = M(3,3) + ZA(iAtom) * (RTx**2+RTy**2)
*
         M(1,2) = M(1,2) + ZA(iAtom) * (    -RTx*RTy)
         M(1,3) = M(1,3) + ZA(iAtom) * (    -RTx*RTz)
         M(2,1) = M(2,1) + ZA(iAtom) * (    -RTy*RTx)

         M(2,3) = M(2,3) + ZA(iAtom) * (    -RTy*RTz)
         M(3,1) = M(3,1) + ZA(iAtom) * (    -RTz*RTx)
         M(3,2) = M(3,2) + ZA(iAtom) * (    -RTz*RTy)
      End Do
*
*     Remove noise
*
      Do i = 1, 3
         Do j = 1, 3
            If (abs(M(i,j)).lt.1.0D-14) M(i,j)=Zero
         End Do
      End Do
C     Call RecPrt('Compute_M: M',' ',M,3,3)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(Z_Tot)
      End
