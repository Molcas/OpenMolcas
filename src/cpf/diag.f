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
      SUBROUTINE DIAG_CPF(ICASE,JSY,HDIAG,FC,FIJ,FJI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JSY(*),HDIAG(*),FC(*),FIJ(*),FJI(*)
      DIMENSION ICASE(*)
      CALL QENTER('DIAG')
      CALL IIJJ(ICASE,JSY,HDIAG,FC,FIJ,FJI)
      CALL IJIJ(JSY,HDIAG,FJI)
      CALL QEXIT('DIAG')
      RETURN
      END
