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
      Subroutine WelMem(nHer,MemWel,la,lb,lr)
*
      k = la+lb
      jsum = 1
      Do 10 i = 1, k
         jsum = jsum + 3**i
 10   Continue
      nHer=1
      MemWel = jsum +Max((k+1)*(k/2+1)*(k/4+1)+1,
     &                   9+3**k,
     &                   5)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
