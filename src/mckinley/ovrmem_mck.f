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
      Subroutine OvrMem_mck(nHer,MmMltP,la,lb,lr)
*
      nElem(i) = (i+1)*(i+2)/2

*
      lr=0
      nHer=(la+lb+lr+3)/2
      MmMltP = 3*nHer*(la+2) +
     &         3*nHer*(lb+2) +
     &         3*nHer*(lr+2) +
     &         3*(la+2)*(lb+2)*(lr+1)+2 +
     &         nelem(la)*nelem(lb)*2
*
      Return
      End
