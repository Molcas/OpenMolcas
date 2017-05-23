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
      SUBROUTINE CHO_CNVTIM(TIM,IHR,IMN,SEC)
C
C     Purpose: convert TIM to hours/minutes/seconds
C
      IMPLICIT NONE
      REAL*8  TIM
      INTEGER IHR, IMN
      REAL*8  XHR, XMN, SEC

      XHR = TIM/3.6D3
      IHR = INT(XHR)

      XMN = (TIM - 3.6D3*DBLE(IHR))/6.0D1
      IMN = INT(XMN)

      SEC = TIM - 3.6D3*DBLE(IHR) - 6.0D1*DBLE(IMN)

      END
