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
      Subroutine ChCar(iChCar,iGen,nGen)
      Implicit Integer (a-z)
      Integer iChCar(3), iGen(nGen)
*
*     Generate characteristics for x, y, and z.
*
      Do iCar = 1, 3
         iChCar(iCar) = 0
         iComp = 2**(iCar-1)
         Do i = 1, nGen
            If (iAnd(iGen(i),iComp).eq.iComp) iChCar(iCar) = iComp
         End Do
      End Do
*
      Return
      End
