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
      Subroutine MltMem(nHer,MemMlt,la,lb,lr)
*
#include "rmat_option.fh"
*
      nElem(i) = (i+1)*(i+2)/2
*
      nHer=(la+lb+lr+2)/2
      nComp = nElem(lr)
      MemMlt = 3*nHer*(la+1) +
     &         3*nHer*(lb+1) +
     &         3*nHer*(lr+1) +
     &         3*(la+1)*(lb+1)*(lr+1) +
     &         nElem(la)*nElem(lb)*nComp
      If (RMat_type_integrals) Then
         MemMlt = MemMlt + la+lb+lr+1
      End If
*
      Return
      End
