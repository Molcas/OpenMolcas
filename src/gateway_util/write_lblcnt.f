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
      Subroutine Write_LblCnt(LuWr,LblCnt,x)
      Implicit Real*8 (A-H,O-Z)
      Character*(*) LblCnt
      Real*8 x(3)
      Write (LuWr,'(1X,A,1X,3F20.10)') LblCnt,x
      Return
      End
