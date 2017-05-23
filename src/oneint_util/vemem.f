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
      Subroutine VeMem(nHer,MemVe,la,lb,lr)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      nHer=(la+(lb+1)+2)/2
      MemVe = 3*nHer*(la+1)
     &      + 3*nHer*(lb+2)
     &      + 3*nHer
     &      + 3*(la+1)*(lb+2)
     &      + 3*(la+1)*(lb+1) + 1
     &      + 3*nElem(la)*nElem(lb)*nELem(lr)
*
      Return
      End
