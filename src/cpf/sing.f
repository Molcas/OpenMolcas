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
C
      SUBROUTINE SING(IWHY)
      IMPLICIT REAL*8 (A-H,O-Z)
11    FORMAT(54H MATRIX WITH ZERO ROW IN DECOMPOSE.                   )
12    FORMAT(54H SINGULAR MATRIX IN DECOMPOSE. ZERO DIVIDE IN SOLVE.  )
13    FORMAT(54H NO CONVERGENCE IN IMPROVE. MATRIX IS NEARLY SINGULAR.)
      NOUT=6
C**** NOUTE=STANDARD OUTPUT UNIT
      GOTO (1,2,3),IWHY
1     WRITE(NOUT,11)
      CALL XFLUSH(6)
      GOTO 10
2     WRITE(NOUT,12)
      CALL XFLUSH(6)
      GOTO 10
3     WRITE(NOUT,13)
      CALL XFLUSH(6)
10    RETURN
      END
