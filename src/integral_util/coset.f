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
      Subroutine CoSet(iCoSet,nCoSet,iChAtom,iOper,nSym)
      Implicit Real*8 (A-H,O-Z)
      Integer iCoSet(0:7), iOper(0:nSym-1)
      Logical Same
*
*     Find the coset representatives
*
      iCoSet(0) = 0      ! Put in the unit operator
      nCoSet = 1
      Do iIrrep = 1, nSym-1
         itest=iAnd(iChAtom,iOper(iIrrep))
         Same = .False.
         Do jCoSet = 0, nCoSet-1
            jTest = iAnd(iChAtom,iCoSet(jCoSet))
            Same = Same .or. jTest.eq.iTest
         End Do
         If (.Not.Same) Then
            nCoSet = nCoSet + 1
            iCoSet(nCoSet-1) = iOper(iIrrep)
         End If
      End Do
*
      Return
      End
