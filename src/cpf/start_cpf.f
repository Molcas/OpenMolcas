************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE START_CPF(C,NCONF,IREF0)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*)
      DO 10 I=1,NCONF
      C(I)=0.0D0
10    CONTINUE
      C(IREF0)=1.0D0
      RETURN
      END
