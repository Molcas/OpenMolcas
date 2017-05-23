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
      SUBROUTINE ZBLTP(ISMOST,MAXSYM,IDC,ICBLTP,IMMLST)
*
* Generate vector ICBLTP giving type of each block
*
*
* ICBLTP gives type of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
*. Input
      DIMENSION ISMOST(*),IMMLST(*)
*. Output
      DIMENSION ICBLTP(*)
*. Changed to simplify structure for IDC .le. 2
      IF(IDC.LE.2) THEN
*. No spatial degeneracy
        DO IASYM = 1, MAXSYM
          IBSYM = ISMOST(IASYM)
          IF(IDC.EQ.2.AND.IBSYM.GT.IASYM) THEN
*.Symmetry block excluded
            ICBLTP(IASYM) = 0
          ELSE IF((IDC.EQ.2.AND.IASYM.GT.IBSYM).OR.IDC.EQ.1) THEN
*.Complete symmetry block included
            ICBLTP(IASYM) = 1
          ELSE
*.Lower half  symmetry block included
            ICBLTP(IASYM) = 2
          END IF
        END DO
      ELSE
*. Also spatial degeneracy
      DO 100 IASYM = 1, MAXSYM
*
        IBSYM = ISMOST(IASYM)
        IF(IBSYM .EQ. 0 ) GOTO 100
        IF(((IDC.EQ.2.OR.IDC.EQ.4).AND.(IBSYM.GT.IASYM))
     &                    .OR.
     &       (IDC.EQ.3.AND.IMMLST(IASYM).GT.IASYM)) THEN
*.Symmetry block excluded
          ICBLTP(IASYM) = 0
        ELSE IF((IDC.EQ.2.AND.IASYM.GT.IBSYM)
     &                   .OR.
     &                IDC.EQ.1
     &                   .OR.
     &          (IDC.EQ.3.AND.IASYM.GE.IMMLST(IASYM))) THEN
*.Complete symmetry block included
          ICBLTP(IASYM) = 1
        ELSE
*.Lower half  symmetry block included
          ICBLTP(IASYM) = 2
        END IF
  100 CONTINUE
      END IF
*     ^ End of IDC switch
*
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Block type of symmetry blocks '
         CALL IWRTMA(ICBLTP,1,MAXSYM,1,MAXSYM)
      END IF
*
      RETURN
      END
