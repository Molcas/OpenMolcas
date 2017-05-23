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
      FUNCTION IMNMX(IVEC,NDIM,MINMAX)
*
*     Find smallest (MINMAX=1) or largest (MINMAX=2)
*     absolute value of elements in integer vector IVEC
*
      DIMENSION IVEC(*)
*
      IX = 0
      IF(NDIM.GT.0) THEN
        IX = -1
        IF(MINMAX.EQ.1) THEN
          IX=ABS(IVEC(1))
          DO I=2,NDIM
            IX=MIN(IX,ABS(IVEC(I)))
          END DO
        END IF
*
        IF(MINMAX.EQ.2) THEN
          IX=ABS(IVEC(1))
          DO I=2,NDIM
            IX=MAX(IX,ABS(IVEC(I)))
          END DO
        END IF
*
      ELSE IF(NDIM.EQ.0) THEN
*. No components : set to zero and write a warning
        IX = 0
        WRITE(6,*) ' Min/Max taken zero length vector set to zero'
      END IF
*
      IMNMX = IX
*
      RETURN
      END
