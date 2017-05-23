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
* Copyright (C) 1990, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE DEGVEC(VEC,NDIM,NDGVL,IDEG)
*
* A vector VEC is given with elements in ascending order
* group elements in degenerate pairs
*
*=======
* Input
*=======
* VEC : input vector
* NDIM : Number of elements in vec
*
*========
* Output
*========
* NDGVL : Number of degenerate values
* IDEG(I) : Number of elements in VEC with degenerate value I
*
* Jeppe Olsen , April 1990
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*.Input
      DIMENSION VEC(*)
*.Output
      DIMENSION IDEG(*)
*.Threshold for defining degenerency
      THRES = 1.0D-8
C?      write(6,*) ' Input vector to DEGVEC '
C?      call wrtmat(VEC,1,NDIM,1,NDIM)
      XDGVL = VEC(1)
      NDEG = 1
      NDGVL = 0
      DO 100 I = 2, NDIM
        IF(ABS(VEC(I)-XDGVL).LE.THRES) THEN
          NDEG = NDEG + 1
        ELSE
          NDGVL = NDGVL + 1
          IDEG(NDGVL) = NDEG
          XDGVL = VEC(I)
          NDEG = 1
        END IF
  100 CONTINUE
*. Last group
      NDGVL = NDGVL + 1
      IDEG(NDGVL) = NDEG
*
      NTEST = 0
      IF(NTEST .GT. 0 ) THEN
        WRITE(6,*) ' Output from DEGVEC '
        WRITE(6,*) ' ================== '
        WRITE(6,*)
        WRITE(6,*) ' Number of degenerate values ' ,NDGVL
        WRITE(6,*) ' Degenerencies of each value '
        CALL IWRTMA(IDEG,1,NDGVL,1,NDGVL)
      END IF
*
      RETURN
      END
