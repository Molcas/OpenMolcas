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
* Copyright (C) 1990, Markus P. Fuelscher                              *
************************************************************************
      SUBROUTINE UPKVEC(NITEM,ICVEC,CVEC)
C***********************************************************************
C
C     PURPOSE:
C     DECODE THE CI-VECTOR BY CHANGING THE NUMBER REPRESENTATION FROM
C     INTEGER*4 ==> REAL*8
C
C**** M.P. FUELSCHER, UNIVERSITY OF LUND, SWEDEN, NOV. 1990 ************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CVEC(NITEM),ICVEC(NITEM)
      PARAMETER (SCALE=1.0D0/2147483647.0D0)
C
      DO 10 ITEM=1,NITEM
        CVEC(ITEM)=SCALE*DBLE(ICVEC(ITEM))
10    CONTINUE
C
      RETURN
      END
