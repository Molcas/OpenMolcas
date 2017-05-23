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
      Subroutine MemRys_g(iSD4,nSD,nRys,MemPrm)
      Implicit Real*8 (a-h,o-z)
      Integer iSD4(0:nSD,4), iAnga(4)
*
      iAnga(1) = iSD4(1,1)
      iAnga(2) = iSD4(1,2)
      iAnga(3) = iSD4(1,3)
      iAnga(4) = iSD4(1,4)
      Call MemRg1(iAnga,nRys,MemPrm)
*
      Return
      End
