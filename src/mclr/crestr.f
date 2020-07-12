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
      SUBROUTINE CRESTR(STRING,NSTINI,NSTINO,NEL,NORB,
     &                  Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,
     &                  ISTMPL,ISTMPO,LROW,I1TYP,IPRNT)
*
* A set of strings containing NEL electrons are given
* set up all possible ways of adding an electron to this set of strings
*
*========
* Input :
*========
*
*
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* Z      : Lexical ordering matrix for output strings containing
*          NEL + 1 electrons
* NEWORD : Reordering array for N+1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL+1 strings
*
* LROW   : Length of tables for each string, negative number
*          indicates compact form
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the resulting string has carries a negative sign
*          then the string number is negative
*          TTO(I,ISTRIN) = LSTRIN.NE.0 indicates that I added to
*                          ISTRIN  gives LSTRIN . If LSTRIN is
*                          gt. NSTINO then
*                          A+(I)!ISTRIN>=-!LSTRIN-NSTINO> .
*. If lrow .le. 0 TI and TTO are constructed in sparse form using
*. ISTMPO and ISTMPL as pointers
* ISTST
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER  STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
C     DIMENSION TI(NORB,NSTINI),TTO(NORB,NSTINI)
      DIMENSION TI(*), TTO(*)
      DIMENSION ISTMPO(*),ISTMPL(*)
*.Scratch
      DIMENSION STRIN2(500)
*
* LCR NOT DECLARED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*  I SET I TO ZERO
      LCR=0
      NTEST0 =  000
      NTEST = MAX(IPRNT,NTEST0)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' CRESTR speaking '
        WRITE(6,*)  ' =============== '
      END IF
      LUOUT = 6
*
      IOFF=0     ! dummy initialize
      IPLACE=-1  ! dummy initialize
      DO 1000 ISTRIN = 1,NSTINI
*. Offset for current creations from this string
       IF(ISTRIN.EQ.1) THEN
         IOFF = 1
       ELSE
         IF(LROW.GT.0) THEN
            IOFF = (ISTRIN-1)*LROW + 1
         ELSE
            IOFF = IOFF + LCR
         END IF
       END IF
       LCR = 0
        DO 100 IORB = 1, NORB
           IF(NEL.EQ.0) THEN
             IPLACE = 1
             GOTO 11
           ELSE IF ( NEL .NE. 0 ) THEN
            DO 10 IEL = 1, NEL
              IF(IEL.EQ.1.AND.STRING(1,ISTRIN).GT.IORB) THEN
                IPLACE = 1
                GOTO 11
              ELSE IF( (IEL.EQ.NEL.AND.IORB.GT.STRING(IEL,ISTRIN)) .OR.
     &                 (IEL.LT.NEL.AND.IORB.GT.STRING(IEL,ISTRIN).AND.
     &                  IORB.LT.STRING(IEL+1,ISTRIN)) ) THEN
                IPLACE = IEL+1
                GOTO 11
              ELSE IF(STRING(IEL,ISTRIN).EQ.IORB) THEN
                IPLACE = 0
                GOTO 11
              END IF
   10       CONTINUE
           END IF
   11     CONTINUE
*
          IF(IPLACE.NE.0) THEN
*. Generate next string
            DO 30 I = 1, IPLACE-1
            STRIN2(I) = STRING(I,ISTRIN)
   30       CONTINUE
            STRIN2(IPLACE) = IORB
            DO 40 I = IPLACE,NEL
            STRIN2(I+1) = STRING(I,ISTRIN)
   40       CONTINUE
*. Is new string allowed ?
            ITYPE = IOCTP2_MCLR(STRIN2,NEL+1,I1TYP)
            IF(ITYPE.NE.0) THEN
              JSTRIN = ISTRNM(STRIN2,NORB,NEL+1,Z,NEWORD,1)
*. Save it !
              IF(LROW.GT.0) THEN
                LCR = IORB
              ELSE
                LCR = LCR + 1
              END IF

              TTO(IOFF-1+LCR ) = JSTRIN
              IIISGN = (-1)**(IPLACE-1)
              IF(LSGSTR.NE.0)
     &        IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
              IF(IIISGN .EQ. -1 ) TTO(IOFF-1+LCR ) = - TTO(IOFF-1+LCR )
              TI(IOFF-1+LCR ) = IORB

C             TTO(IORB,ISTRIN) = JSTRIN
C             IIISGN = (-1)**(IPLACE-1)
C             IF(LSGSTR.NE.0)
C    &        IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
C             IF(IIISGN .EQ. -1 ) TTO(IORB,ISTRIN) = - TTO(IORB,ISTRIN)
C             TI(IORB,ISTRIN) = IORB
            END IF
          END IF
  100   CONTINUE
*
        IF(LROW.LT.0) THEN
          ISTMPO(ISTRIN) = IOFF
          ISTMPL(ISTRIN) = LCR
C          write(6,*) ' ISTRIN ISTMPO ISTMPL '
C          write(6,*) ISTRIN, ISTMPO(ISTRIN), ISTMPL(ISTRIN)
        END IF
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(LUOUT,*) ' Output from CRESTR : '
        WRITE(LUOUT,*) '==================='
*
        IF(LROW.GT.0) THEN
          WRITE(6,*) ' Full map '
          WRITE(6,*)
          WRITE(LUOUT,*) ' Strings with an electron added  '
          DO ISTRIN = 1, NPR
             WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &       'String..',ISTRIN,' New strings.. ',
     &       (TTO((ISTRIN-1)*LROW+I),I = 1,NORB)
          END DO
          DO ISTRIN = 1, NPR
             WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &       'String..',ISTRIN,' orbitals added or removed ' ,
     &       (TI((ISTRIN-1)*LROW+I),I = 1,NORB)
          END DO
        ELSE
          WRITE(6,*) ' Compact map '
          WRITE(6,*)
          WRITE(LUOUT,*) ' Strings with an electron added  '
          DO ISTRIN = 1, NPR
             WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &       'String..',ISTRIN,' New strings.. ',
     &       (TTO(ISTMPO(ISTRIN)-1+I),I = 1,ISTMPL(ISTRIN))
          END DO
          DO ISTRIN = 1, NPR
             WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &       'String..',ISTRIN,' orbitals added or removed ' ,
     &       (TI(ISTMPO(ISTRIN)-1+I),I = 1,ISTMPL(ISTRIN))
          END DO
        END IF
      END IF
*
      RETURN
      END
