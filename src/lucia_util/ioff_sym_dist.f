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
      FUNCTION IOFF_SYM_DIST(ISYM,NGASL,IOFF,MAXVAL,MINVAL)
*
* A ts block of string is given and the individual
* symmetrydisrtributions has been obtained ( for example
* by TS_SYM_PNT)
*
* Obtain offset for symmetrycombination defined by ISYM
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ISYM(*),IOFF(*),MAXVAL(*),MINVAL(*)
* Address in IOFF is
*     1
*     +  (ISM1-MINVAL(1))
*     +  (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     +  (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
*
C     write(6,*) ' Isym and minval '
C     call iwrtma(isym,1,ngasl,1,ngasl)
C     call iwrtma(minval,1,ngasl,1,ngasl)
      I = 1
      IMULT = 1
      DO IGAS = 1, NGASL-1
        I = I + (ISYM(IGAS)-MINVAL(IGAS)) * IMULT
        IMULT = IMULT*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
C       write(6,*) ' igas i imult ',igas,i,imult
      END DO
      IOFF_SYM_DIST=IOFF(I)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from IOFF_SYM_DIST'
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Address and offset ',I,IOFF_SYM_DIST
        WRITE(6,*) ' Symmetry distribution : ', (ISYM(J),J=1,NGASL)
      END IF
*
      RETURN
      END
