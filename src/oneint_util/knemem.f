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
      Subroutine KnEMem(nHer,MemKne,la,lb,lr)
*
#include "rmat_option.fh"
*
      nHer=(la+lb+lr+2)/2
      MemKnE = 3*nHer*(la+2) +
     &         3*nHer*(lb+2) +
     &         3*nHer*(lr-1) +
     &         3*(la+2)*(lb+2)*(lr-1) +
     &         3*(la+1)*(lb+1) + 1 + 1
      If (RMat_type_integrals) Then
         MemKne = MemKne
     &          + la+lb+3
     &          + la+lb+1
     &          + la+lb+1
      End If
*
      Return
      End
