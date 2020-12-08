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
      Subroutine NACInt(xyz,nCent,H12,Bf,lWrite_,Label,dBf,ldB)
      use Slapaf_Info, only: NAC
      Implicit Real*8  (a-h,o-z)
#include "info_slapaf.fh"
#include "real.fh"
#include "nadc.fh"
#include "constants.fh"
      Real*8   Bf(3,nCent), xyz(3,nCent), dBf(3*nCent,3*nCent)
      Logical lWrite_, ldB
      Character*8 Label
*
*
      H12=Zero
      If (lWrite_) Then
         Write (6,'(2A,F18.8,A,F18.8,A)')
     &             Label,' : H12               = ',
     &             H12, ' hartree '
      End If
*
*---- Compute the WDC B-matrix
*
C     Call RecPrt('NAC',' ',NAC,3,nCent)
      Do iCent = 1, nCent
         Fact=DBLE(iDeg(xyz(1,iCent)))
C        Write (6,*) 'Fact=',Fact
         Do iCar = 1, 3
            Bf(iCar,iCent)=Fact*NAC(iCar,iCent)
         End Do
      End Do
C     Call RecPrt('Bf',' ',Bf,3,nCent)
*
*---- Compute the cartesian derivative of the B-Matrix.
*
      If (ldB) Then
         nGrad=3*nCent
         Call FZero(dBf,(3*nCent)**2)
*        Call RecPrt('dBf','(9F9.1)',dBf,3*nCent,3*nCent)
*
      End If
*
      Return
      End
