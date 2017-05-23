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
C     Integer Function Set_CHO_ADRVEC(ii)
      Function Set_CHO_ADRVEC(ii)
#include "cholesky.fh"
      Integer Set_CHO_ADRVEC
*
      Set_CHO_ADRVEC=0
      If (ii.lt.0) Then
         Set_CHO_ADRVEC=CHO_ADRVEC
      Else If (ii.eq.1.or.ii.eq.2) Then
         CHO_ADRVEC=ii
         Set_CHO_ADRVEC=CHO_ADRVEC
      Else
         Call WarningMessage(2,'Set_CHO_ADRVEC: Illegal option')
         Write (6,*) 'ii=',ii
         Call Abend()
      End if
*
      Return
      End
