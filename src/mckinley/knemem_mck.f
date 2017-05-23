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
      Subroutine Knemem_mck(nHer,MmKnEG,la,lb,lr)
*
      nElem(i) = (i+1)*(i+2)/2
      nHer=((la+1)+(lb+1)+2+1)/2
      MmKnEG =  3*nHer*(la+3) +
     &         3*nHer*(lb+3) +
     &         3*nHer +
     &         3*(la+3)*(lb+3) +
     &         3*(la+2)*(lb+2) + 1 + 1 +
     &         nElem(la)*nElem(lb)*3
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
