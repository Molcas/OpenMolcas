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
      SUBROUTINE ANNSTR_GAS( STRING, NSTINI, NSTINO,    NEL,   NORB,
     &                       IORBOF,      Z, NEWORD, LSGSTR, ISGSTI,
     &                       ISGSTO,     TI,    TTO,  NACOB,    IEC,
     &                         LDIM,  IPRNT)
*
* A group of strings containing NEL electrons is given
* Set up all possible ways of removing an electron
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* IORBOF : Number of first orbital
* Z      : Lexical ordering matrix for output strings containing
*          NEL - 1 electrons
* NEWORD : Reordering array for N-1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL-1 strings
* IEC    : = 1 Extended map, dimension equals number of orbs
* IEC    : = 2 Compact  map, dimension equals number of elecs
* LDIM   : Row dimension ( see IEC)
*
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the string have a negative sign
*          then the phase equals - 1
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(LDIM,NSTINI),TTO(LDIM,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTEST0 =  1
      NTEST = MAX(IPRNT,NTEST0)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' ANNSTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
*
      DO 1000 ISTRIN = 1,NSTINI
        DO 100 IEL = 1, NEL
*. String with electron removed
          DO JEL = 1, IEL-1
           STRIN2(JEL) = STRING(JEL,ISTRIN)
          END DO
          DO JEL = IEL+1, NEL
           STRIN2(JEL-1) = STRING(JEL,ISTRIN)
          END DO
          JSTRIN = ISTRNM(STRIN2,NACOB,NEL-1,Z,NEWORD,1)
C?        WRITE(6,*) ' anni-string and number '
C?        CALL IWRTMA(STRIN2,1,NEL-1,1,NEL-1)
C?        WRITE(6,*) ' JSTRIN = ', JSTRIN
*
          IORBABS = STRING(IEL,ISTRIN)
          IORB = STRING(IEL,ISTRIN)-IORBOF+1
          IF(IEC.EQ.1) THEN
            IROW = IORB
          ELSE
            IROW = IEL
          END IF
*
          TI(IROW,ISTRIN ) = -IORBABS
C         TI(IROW,ISTRIN ) = -IORB
          TTO(IROW,ISTRIN) = JSTRIN
*PAM2009          IIISGN = (-1)**(IEL-1)
          IIISGN = 1-2*MOD(IEL+1,2)
          IF(LSGSTR.NE.0)
     &    IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
          IF(IIISGN .EQ. -1 )
     &    TTO(IROW,ISTRIN) = - TTO(IROW,ISTRIN)
  100   CONTINUE
*
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(6,*) ' Output from ANNSTR : '
        WRITE(6,*) '==================='
*
        WRITE(6,*)
        WRITE(6,*) ' Strings with an electron added or removed'
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' New strings.. ',
     &     (TTO(I,ISTRIN),I = 1,LDIM)
        END DO
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' orbitals added or removed ' ,
     &     (TI(I,ISTRIN),I = 1,LDIM)
        END DO
      END IF
*
      RETURN
      END
*
*
