!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      INTEGER FUNCTION ISYMS1(STRING,NEL)
!
! Symmmetry of string, D2H version
!
      use symmetry_info, only: SYMPRO => Mul
      use lucia_data, only: ISMFTO
      IMPLICIT None
!. Specific input
      Integer NEL
      INTEGER STRING(*)

      INTEGER ISYM,IEL,NTEST
!
      ISYM = 1
      DO 100 IEL = 1, NEL
        ISYM = SYMPRO(ISYM,ISMFTO(STRING(IEL)))
  100 CONTINUE
!
      ISYMS1 = ISYM
!
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
!
      END FUNCTION ISYMS1
