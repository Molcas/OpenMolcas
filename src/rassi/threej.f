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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
      REAL*8 FUNCTION THREEJ(XJ1,XJ2,XJ3,XM1,XM2,XM3)
      IMPLICIT REAL*8 (A-H,O-Z)
C ThreeJ: REAL*8 Wigner 3-j coefficients. From a modification of
C Racah''s formula for Clebsch-Gordan coeffs.

      THREEJ=DCLEBS(XJ1,XJ2,XJ3,XM1,XM2,-XM3)
      IF(THREEJ.EQ.0.0D0) RETURN
      I=NINT(XJ1-XJ2-XM3)
      IF(I.NE.(I/2)*2) THREEJ=-THREEJ
      THREEJ=THREEJ/SQRT(2.0D0*XJ3+1.0D0)

      RETURN
      END
