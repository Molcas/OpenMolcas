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
      SUBROUTINE MKDAW_RASSI(NLEV,NVERT,IDRT,IDOWN,IDAW,LTV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4),IDRT(NVERT,5)
      DIMENSION LTV(-1:NLEV)

      CALL MKLTV(NVERT,NLEV,IDRT,LTV)

      CALL MKDAW(NVERT,IDOWN,IDAW)

      END
