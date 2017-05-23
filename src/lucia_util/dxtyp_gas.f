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
      SUBROUTINE DXTYP_GAS(   NDXTP,     ITP,     JTP,     KTP,     LTP,
     &                        NOBTP,      IL,      IR)
*
* Obtain types of I,J,K,l so
* <L!a+I a+K a L a J!R> is nonvanishing
* only combinations with type(I) .ge. type(K) and type(J).ge.type(L)
* are included
*
      INTEGER IL(NOBTP),IR(NOBTP)
      INTEGER ITP(*),JTP(*),KTP(*),LTP(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DXTYP_GAS in action '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' Occupation of left string '
        CALL IWRTMA(IL,1,NOBTP,1,NOBTP)
        WRITE(6,*) ' Occupation of right string '
        CALL IWRTMA(IR,1,NOBTP,1,NOBTP)
      END IF
*
*. Number of differing occupations
      NANNI = 0
      NCREA = 0
      NDIFT = 0
*
      ICREA1 = 0
      ICREA2 = 0
      IANNI1 = 0
      IANNI2 = 0
      DO IOBTP = 1, NOBTP
        NDIFT = NDIFT + ABS(IL(IOBTP)-IR(IOBTP))
        NDIF = IL(IOBTP)-IR(IOBTP)
        IF(NDIF.EQ.2) THEN
*. two electrons of type IOBTP must be created
          ICREA1 = IOBTP
          ICREA2 = IOBTP
          NCREA = NCREA + 2
        ELSE IF (NDIF .EQ. -2 ) THEN
*. Two electrons of type IOBTP must be annihilated
          IANNI1 = IOBTP
          IANNI2 = IOBTP
          NANNI = NANNI + 2
        ELSE IF (NDIF.EQ.1) THEN
*. one electron of type IOBTP must be created
          IF(NCREA.EQ.0) THEN
            ICREA1 = IOBTP
          ELSE
            ICREA2 = IOBTP
          END IF
          NCREA = NCREA + 1
        ELSE IF (NDIF.EQ.-1) THEN
* One electron of type IOBTP must be annihilated
          IF(NANNI.EQ.0) THEN
            IANNI1 = IOBTP
          ELSE
            IANNI2 = IOBTP
          END IF
          NANNI = NANNI + 1
        END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NCREA, NANNI ', NCREA, NANNI
        WRITE(6,*)  ' ICREA1, IANNI1 ', ICREA1,IANNI1
        WRITE(6,*)  ' ICREA2, IANNI2 ', ICREA2,IANNI2
      END IF
*
      NDXTP = 0
      IF(NDIFT.GT.4) THEN
        NDXTP = 0
      ELSE
      IF(NCREA.EQ.0.AND.NANNI.EQ.0) THEN
*. strings identical, include diagonal excitions  itp = jtp, ktp=ltp
        DO IJTP = 1, NOBTP
          IF(IR(IJTP).GE.1) THEN
            DO KLTP = 1, IJTP
              IF((IJTP.NE.KLTP.AND.IR(KLTP).GE.1).OR.
     &           (IJTP.EQ.KLTP.AND.IR(KLTP).GE.2)) THEN
                 NDXTP = NDXTP + 1
                 ITP(NDXTP) = IJTP
                 JTP(NDXTP) = IJTP
                 KTP(NDXTP) = KLTP
                 LTP(NDXTP) = KLTP
              END IF
            END DO
          END IF
        END DO
*. Strings differ by single excitation
      ELSE IF( NCREA.EQ.1.AND.NANNI.EQ.1) THEN
*. diagonal excitation plus creation in ICREA1 and annihilation in IANNI1
        DO IDIA = 1, NOBTP
          IF((IDIA.NE.IANNI1.AND.IR(IDIA).GE.1).OR.
     &       (IDIA.EQ.IANNI1.AND.IR(IDIA).GE.2)) THEN
             NDXTP = NDXTP + 1
             ITP(NDXTP) = MAX(ICREA1,IDIA)
             KTP(NDXTP) = MIN(ICREA1,IDIA)
             JTP(NDXTP) = MAX(IANNI1,IDIA)
             LTP(NDXTP) = MIN(IANNI1,IDIA)
          END IF
        END DO
      ELSE IF(NCREA.EQ.2.AND.NANNI.EQ.2) THEN
*. Strings differ by double excitation
        NDXTP = 1
        ITP(1) = ICREA2
        KTP(1) = ICREA1
        JTP(1) = IANNI2
        LTP(1) = IANNI1
      END IF
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,I4)')
     &  ' Number of connecting double excitations ', NDXTP
        IF(NDXTP.NE.0) THEN
          WRITE(6,*) '  ITYP KTYP LTYP JTYP '
          WRITE(6,*) '  ===================='
          DO  IDX = 1,NDXTP
            WRITE(6,'(1H ,4I5)')ITP(IDX),KTP(IDX),LTP(IDX),JTP(IDX)
          END DO
        END IF
      END IF
*
      RETURN
      END
