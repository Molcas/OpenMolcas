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
      Subroutine PPMem(nHer,MemPP,la,lb,lr)
*
      nElem(i) = (i+1)*(i+2)/2
*
      nHer=0
      MemPP=0
      intmax=Max(nElem(la),nElem(lb))
      intmax=intmax**2
      MemPP=MemPP+3*intmax
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
