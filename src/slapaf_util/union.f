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
      SubRoutine Union(iU,nU,iV,nV,iR,iM,nM)
      Implicit Real*8 (A-H,O-Z)
      Integer   iU(nU),  iV(nV), iM(8)
      Logical RinT_
*
*        M is formed as U union RU
*
      Call iCopy(nU,iU,1,iM,1) ! copy the first elements
      nM = nU
      Do i = 1, nV
         iRV = iEor(iR,iV(i))
         If (.Not.RinT_(iM,nM,iRV)) Then
            nM = nM + 1
            iM(nM) = iRV
         End If
      End Do
*
      Return
      End
