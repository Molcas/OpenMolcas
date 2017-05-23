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
      SUBROUTINE ANNSTR(STRING,NSTINI,NSTINO,NEL,NORB,
     &                  Z,NEWORD,LROW,LSGSTR,ISGSTI,ISGSTO,TI,TTO,
     &                  I1TYP,IPRNT)
*
* A set of strings containing NEL electrons are given
* set up all possible ways of annihilating an electron from
* this set of string
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of NEL  strings
* NSTINO : Number of NEL-1 strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* Z      : Lexical ordering matrix for output strings containing
*          NEL - 1 electrons
* NEWORD : Reordering array for N-1 strings
* LROW   : Number of rows in output vector
*          = NEL : Information is written in truncated form
*                  row  corresponds to place of electron
*          = NORB: Information is written in expanded form,
*                  row corresponds to full orbital number
* LSGSTI : NE.0 => use ISGSTI,ISGSTO to allow for sign of string
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL-1 strings
*
*=========
* Output :
*=========
*TI      : Array giving minus orbital annihilated
*          TI(I,ISTRIN).LT.0 : Orbital TI(I,ISTRIN) can be
*          annihilated from string ISTRIN
*TTO     : Resulting NEL - 1 strings
*          if the resulting string has carries a negative sign
*          then the string number is shifted with a NSTINO
*          TTO(I,ISTRIN) = .LE.0 indicates that orbital I cannot be
*                             added to istrin
*          TTO(I,ISTRIN) = LSTRIN.NE.0 indicates that I added to
*                          ISTRIN  gives LSTRIN . If LSTRIN is
*                          nehative
*                          A+(I)!ISTRIN>=-!-LSTRIN> .
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER  STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL-1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(LROW,NSTINI),TTO(LROW,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTEST0 =   0
      NTEST = MAX(NTEST0,IPRNT)
      IF( NTEST .GE.  20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' ANNSTR speaking '
        WRITE(6,*)  ' =============== '
      END IF
      LUOUT = 6
*. Expanded or truncated form
      IF(LROW.EQ.NEL.AND.NEL.NE.NORB)THEN
        IEXPN = 0
      ELSE
        IEXPN = 1
      END IF
*. Loop over input strings
      DO 1000 ISTRIN = 1,NSTINI
*. loop over electrons to be removed
        DO 100 IEL = 1,NEL
          IF(IEXPN.EQ.0) THEN
            IPLACE = IEL
          ELSE
            IPLACE = STRING(IEL,ISTRIN)
          END IF
          DO 30 I = 1, IEL-1
   30     STRIN2(I) = STRING(I,ISTRIN)
          DO 40 I = IEL+1,NEL
   40     STRIN2(I-1) = STRING(I,ISTRIN)
*. Is new string allowed ?
          ITYPE = IOCTP2_MCLR(STRIN2,NEL-1,I1TYP)
          IF(ITYPE.NE.0) THEN
*                    ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
            JSTRIN = ISTRNM(STRIN2,NORB,NEL-1,Z,NEWORD,1)
            TTO(IPLACE,ISTRIN) = JSTRIN
            TI(IPLACE,ISTRIN) = - STRING(IEL,ISTRIN)
            IIISGN = (-1)**(IEL-1)
            IF(LSGSTR.GT.0)
     &      IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
            IF(IIISGN .EQ. -1 ) TTO(IPLACE,ISTRIN) = -TTO(IPLACE,ISTRIN)
          END IF
  100   CONTINUE
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(LUOUT,*) ' Output from ANNSTR : '
        WRITE(LUOUT,*) '==================='
         IF(IEXPN.EQ.0) THEN
           WRITE(LUOUT,*) ' Strings with an electron removed '
         ELSE
           WRITE(LUOUT,*) ' Combined N+1/N-1 string array '
         END IF
         DO 1235 ISTRIN = 1, NPR
            WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &      'String..',ISTRIN,' New strings.. ',
     &      (TTO(I,ISTRIN),I = 1,LROW)
 1235    CONTINUE
*
         WRITE(6,*) ' orbitals removed '
         DO 1236 ISTRIN = 1, NPR
            WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &      'String..',ISTRIN,' orbitals annihilated.. ',
     &      (TI(I,ISTRIN),I = 1,LROW)
 1236    CONTINUE

      END IF
*
      RETURN
      END
