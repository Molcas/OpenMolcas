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
      SUBROUTINE IVCSUM(      IA,      IB,      IC,   IFACB,   IFACC,
     &                      NDIM)
*
* Add two (scaled) integer vectors
*
*        IA(*) = IFACB*IB(*) + IFACC*IC(*)
*
      DIMENSION IA(*),IB(*),IC(*)
*
      DO 100 I = 1, NDIM
        IA(I) = IFACB * IB(I) + IFACC * IC(I)
  100 CONTINUE
*
      RETURN
      END
