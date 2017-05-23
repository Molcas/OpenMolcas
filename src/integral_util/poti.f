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
      Subroutine poti(k,ipot3)
      Implicit Real*8(A-H,O-Z)
      Integer ipot3(0:15)
*
      isum=1
      ipot3(0)=1
      Do 10 i=1,k
         ipot3(i)=3**i
         isum=isum+ipot3(i)
10    Continue
      ipot3(k+1)=isum
*
      Return
      End
