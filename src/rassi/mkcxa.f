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
      SUBROUTINE MKCXA(NSYM,NOSH,NCXA,TRA,CXA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TRA(NCXA),CXA(NCXA)
      DIMENSION NOSH(NSYM)
      ISTA=1
      DO I=1,NSYM
        NDIMEN=NOSH(I)
        IF(NDIMEN.GT.0) THEN
          CALL MKCXAL(NDIMEN,TRA(ISTA),CXA(ISTA))
          ISTA=ISTA+NDIMEN**2
        END IF
      END DO
      RETURN
      END
