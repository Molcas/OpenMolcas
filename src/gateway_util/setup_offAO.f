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
      Subroutine Setup_OffAO()
      use Basis_Info, only: nCnttp, Shells, dbsc
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
*
*     For some reason we need a counter which for a given shell index,
*     kSh, counts how many angular functions have been in the
*     proceeding shells. This is stored in kOffAO(kSh).
*     Additionally, lOffAO gives the total number of angular functions
*     in a given dbsc.
*
      Do iCnttp = 1, nCnttp
         lComp = 0
         Do lSh = 0, nVal_Shells(iCnttp)-1
            kSh = lSh + ipVal(iCnttp)
            If (Shells(kSh)%Prjct ) Then
               kComp = 2*lSh + 1
            Else
               kComp = (lSh+1)*(lSh+2)/2
            End If
            Shells(kSh)%kOffAO = lComp
            If (Shells(kSh)%nBasis_C.ne.0 .and.
     &          Shells(kSh)%nExp    .ne.0 ) lComp = lComp + kComp
         End Do
         dbsc(iCnttp)%lOffAO = lComp
      End Do
*
      Return
      End
