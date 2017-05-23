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
      Subroutine Gen_iSD4(iS, jS, kS, lS,iSD,nSD,iSD4)
      Implicit Real*8 (a-h,o-z)
      Integer iSD(0:nSD,1024), iSD4(0:nSD,4), jQuad(4)
*
      jQuad(1)=iS
      jQuad(2)=jS
      jQuad(3)=kS
      jQuad(4)=lS
      Do i = 1, 4
         iSkal=jQuad(i)
         Do j = 0, nSD
            iSD4(j,i)=iSD(j,iSkal)
         End Do
      End Do
*
      Return
      End
