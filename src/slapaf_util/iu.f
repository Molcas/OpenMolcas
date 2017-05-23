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
      Integer Function iU(iStab,nStab)
      Integer iStab(0:nStab-1)
      iU=0
      Do i = 0, nStab-1
         iU = iOr(iU,2**iStab(i))
      End Do
      Return
      End
