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
      SUBROUTINE ZOOS(ISMOST,IBLTP,MAXSYM,IOCOC,NSSOA,NSSOB,
     &                NOCTPA,NOCTPB,IDC,IOOS,NOOS,NCOMB,IXPND)
*
* Generate offsets for CI vector for RAS CI expansion of given symmetry
* Combination type is defined by IDC
* Total number of combinations NCOMB is also obtained
*
* Symmetry is defined through ISMOST
*
* ICBLTP gives typo of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
* If IXPND .ne. 0 , the diagonal blocks are always choosen expanded
*
* ========
*  Output
* ========
*
* IOOS(IOCA,IOCB,ISYM) : Start of block with alpha strings of
*                        symmetry ISYM and type IOCA, and
*                        betastrings of type IOCB
* NOOS(IOCA,IOCB,ISYM) : Number of combinations
* The ordering used for the CI vector is
*
*    SYMMETRY  FOR ALPHA STRINGS..(GIVES SYMMETRY OF BETA STRING )
*         OCCUPATION TYPE  FOR ALPHA STRING
*            OCCUPATION TYPE FOR    BETA STRING
*                BETA STRING ( COLUMN INDEX)
*                ALPHA STRINGS ( ROW INDEX )
*    END OF LOOPS
*
*
*. Input
      DIMENSION IOCOC(NOCTPA,NOCTPB),ISMOST(*)
      DIMENSION NSSOA(MAXSYM,NOCTPA),NSSOB(MAXSYM,NOCTPB)
      DIMENSION IBLTP(*)
*. output
      DIMENSION IOOS(NOCTPA,NOCTPB,MAXSYM)
      DIMENSION NOOS(NOCTPA,NOCTPB,MAXSYM)
*
      CALL ISETVC(IOOS,0,NOCTPA*NOCTPB*MAXSYM)
      CALL ISETVC(NOOS,0,NOCTPA*NOCTPB*MAXSYM)
C?    CALL ISETVC(ICBLTP,0,MAXSYM)
      NCOMB = 0
      DO 100 IASYM = 1, MAXSYM
*
        IBSYM = ISMOST(IASYM)
        IF(IBSYM .EQ. 0 ) GOTO 100
*. Allowed combination symmetry block ?
        IF(IDC.NE.1.AND.IBLTP(IASYM).EQ.0) GOTO 100
*. Allowed occupation combinations
        DO  95 IAOCC = 1,NOCTPA
          IF(IBLTP(IASYM).EQ.1) THEN
            MXBOCC = NOCTPB
            IREST1 = 0
          ELSE
            MXBOCC = IAOCC
            IREST1 = 1
          END IF
          DO 90 IBOCC = 1, MXBOCC
*.Is this block allowed
            IF(IOCOC(IAOCC,IBOCC).EQ.1) THEN
              IOOS(IAOCC,IBOCC,IASYM) = NCOMB+1
              IF(IXPND.EQ.0 .AND. IREST1.EQ.1 .AND. IAOCC.EQ.IBOCC)THEN
                NCOMB = NCOMB
     &      +   (NSSOA(IASYM,IAOCC)+1)*NSSOB(IBSYM,IBOCC)/2
                NOOS(IAOCC,IBOCC,IASYM) =
     &          (NSSOA(IASYM,IAOCC)+1)*NSSOB(IBSYM,IBOCC)/2
              ELSE
                NCOMB = NCOMB
     &      +   NSSOA(IASYM,IAOCC)*NSSOB(IBSYM,IBOCC)
                NOOS(IAOCC,IBOCC,IASYM) =
     &          NSSOA(IASYM,IAOCC)*NSSOB(IBSYM,IBOCC)
              END IF
            END IF
!            write(6,*) ' NOOS(IA,IB,ISM) ',NOOS(IAOCC,IBOCC,IASYM)
   90     CONTINUE
   95   CONTINUE
  100 CONTINUE
*


      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*)
         WRITE(6,*) ' ==============='
         WRITE(6,*) ' ZOOS reporting '
         WRITE(6,*) ' ==============='
         WRITE(6,*)
         WRITE(6,*) ' NSSOA, NSSOB ( input ) '
         CALL IWRTMA(NSSOA,MAXSYM,NOCTPA,MAXSYM,NOCTPA)
         CALL IWRTMA(NSSOB,MAXSYM,NOCTPB,MAXSYM,NOCTPB)
         WRITE(6,*)
         WRITE(6,*) ' Number of combinations obtained ',NCOMB
         WRITE(6,*) ' Offsets for combination OOS blocks '
         DO 111 IASYM = 1,MAXSYM
           WRITE(6,'(A,I2)') '  Symmetry ',IASYM
           CALL IWRTMA(IOOS(1,1,IASYM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  111    CONTINUE
         WRITE(6,*) ' Number of  combinations per OOS blocks '
         DO 112 IASYM = 1,MAXSYM
           WRITE(6,'(A,I2)') '  Symmetry ',IASYM
           CALL IWRTMA(NOOS(1,1,IASYM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  112    CONTINUE
      END IF
*
      RETURN
      END
