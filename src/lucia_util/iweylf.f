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
      FUNCTION IWEYLF(NOPEN,MULTS)
C
C NUMBER OF CSF'S WITH NOPEN ORBITALS AND TOTAL MULTIPLICITY
C MULTS ACCORDING TO WEYLS FORMULAE
C
C     (2S+1)/(NOPEN+1) * BION(NOPEN+1/0.5NOPEN-S)
C
      IMPLICIT REAL*8           (A-H,O-Z)
C
      NTEST = 00

      IF(NOPEN.EQ.0 .AND. MULTS .EQ. 1 ) THEN
        NCSF = 1
      ELSEIF(MOD(MULTS-1,2) .NE. MOD(NOPEN,2) ) THEN
        NCSF = 0
      ELSEIF(MOD(MULTS-1,2) .EQ. MOD(NOPEN,2) ) THEN
        NCSF = MULTS*IBION_LUCIA(NOPEN+1,(NOPEN+1-MULTS)/2)/(NOPEN+1)
      END IF
C
      IWEYLF = NCSF
C
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,'(A,4I4)')
     &  '  IWEYLF SAYS : NOPEN MULTS NCSF : ', NOPEN,MULTS,NCSF
      END IF
C
      RETURN
      END

*
