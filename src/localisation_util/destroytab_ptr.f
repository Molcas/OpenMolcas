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
      Subroutine DestroyTab_ptr(nAtoms,nOrb2loc,iTab_ptr)
      Implicit Real*8(a-h,o-z)
      Integer iTab_ptr(*)
c
      Do iAt=1,nAtoms
        ip=iTab_ptr(iAt)
        Call Getmem('PA__dum','Free','Real',ip,nOrb2loc**2)
      End Do
c
      Return
      End
