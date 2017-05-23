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
      SUBROUTINE SIGVST(ISGVST,NSMST)
*
* Obtain ISGVST(ISM) : Symmetry of sigma v on string of symmetry ism
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ISGVST(*)
*
      DO 100 ISM = 1, NSMST
C            MLSM(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM(IML,IPARI,ISM,'ST',2)
        MIML = - IML
        CALL MLSM(MIML,IPARI,MISM,'ST',1)
        ISGVST(ISM) = MISM
  100 CONTINUE
*
      NTEST = 1
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' ISGVST array '
        WRITE(6,*) ' ============ '
        CALL IWRTMA(ISGVST,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
