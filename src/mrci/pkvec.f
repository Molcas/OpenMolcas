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
      SUBROUTINE PKVEC(NITEM,CVEC,ICVEC)
C***********************************************************************
C
C     PURPOSE:
C     ENCODE THE CI-VECTOR BY CHANGING THE NUMBER REPRESENTATION FROM
C     REAL*8 ==> INTEGER*4
C
C     NOTE:
C     THE INCOMING DATA CVEC SHOULD NOT BE GREATER THAN 1.0.
C     THE ACCURACY OF THE UNPACKED VALUES IS APPROX. 1.0E-9.
C
C**** M.P. FUELSCHER, UNIVERSITY OF LUND, SWEDEN, NOV. 1990 ************
C
      IMPLICIT REAL*8 (A-H,O-Z)
* NOTE VERY CAREFULLY!! NINT() (and maybe similar) intrinsics
* are BROKEN in CYGWIN GFORTAN!!
* So the intermediate variable tmp is necessary for it to work!
#ifdef _CYGWIN_
      REAL*4 tmp
#endif
      DIMENSION CVEC(NITEM),ICVEC(NITEM)
      PARAMETER (SCALE=2147483647.0D0)
      INTRINSIC NINT
C
      DO 10 ITEM=1,NITEM
#ifdef _CYGWIN_
         tmp=SCALE*CVEC(ITEM)
         ICVEC(ITEM)=NINT(tmp)
#else
        ICVEC(ITEM)=NINT(SCALE*CVEC(ITEM))
#endif
10    CONTINUE
C
      RETURN
      END
