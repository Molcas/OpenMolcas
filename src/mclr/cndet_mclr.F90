!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************
      SUBROUTINE CNDET_MCLR(ICONF,IPDET,NDET,NEL,NORB,NOP,NCL,          &
     &                 IDET,IPRNT)
!
! A configuration ICONF in compressed form and a set of
! prototype determinants ,IPDET, are given .
!
! Construct the corresponding determinants in contracted  form .
!
! JEPPE OLSEN , NOVEMBER 1988
!
      IMPLICIT NONE
      Integer NEL
      Integer ICONF(*   )
!     Integer IPDET(NOP,NDET)
      Integer IPDET(*       )
      Integer NORB,NOP,NCL
      Integer IDET(NEL,*   )
      Integer IPRNT

! local variables
      Integer NTEST,ICL,IBASE,JDET,NDET,IADD,IOP,IADR
!
!
! POSITIVE NUMBER  : ALPHA ORBITAL
! NEGATIVE NUMBER  : BETA  ORBITAL
!
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
      IF( NTEST .GT.200 ) THEN
        IF(NCL .NE. 0 ) THEN
          WRITE(6,*) ' DOUBLE OCCUPIED ORBITALS '
          CALL IWRTMA(ICONF,1,NCL,1,NCL)
        END IF
        IF(NOP .NE. 0 ) THEN
          WRITE(6,*) ' OPEN ORBITALS '
          CALL IWRTMA(ICONF(1+NCL),1,NOP,1,NOP)
        END IF
      END IF
!
!.. 1 DOUBLY OCCUPIED ORBITALS ARE PLACED FIRST
!
      DO 100 ICL = 1, NCL
        IBASE = 2 * (ICL-1)
        DO 90 JDET = 1, NDET
          IDET(IBASE+1,JDET) =  ICONF(ICL)
          IDET(IBASE+2,JDET) = -ICONF(ICL)
90      CONTINUE
100   CONTINUE
!
!..2  SINGLY OCCUPIED ORBITALS
!
      IADD = 2*NCL
      DO 200 JDET = 1, NDET
        DO 190 IOP = 1, NOP
          IADR = (JDET-1)*NOP + IOP
          IF( IPDET(IADR    ) .EQ. 1 ) IDET(IADD+IOP,JDET) =            &
     &    ICONF(NCL +IOP)
          IF( IPDET(IADR    ) .EQ. 0 ) IDET(IADD+IOP,JDET) =            &
     &    - ICONF(NCL +IOP)
190     CONTINUE
200   CONTINUE
!
!..3 OUTPUT
!
      IF( NTEST.GE.200) THEN
       WRITE(6,*) ' CONFIGURATION FROM DETCON '
       CALL IWRTMA(ICONF,1,NORB,1,NORB)
       WRITE(6,* ) ' PROTO TYPE DETERMINANTS '
       IF(NOP*NDET .GT. 0)                                              &
     & CALL IWRTMA(IPDET,NOP,NDET,NOP,NDET)
       IF(NEL*NDET .GT. 0 )                                             &
     & WRITE(6,*) ' CORRESPONDING DETERMINANTS '
       CALL IWRTMA(IDET,NEL,NDET,NEL,NDET)
      END IF
!
!..4  EXIT
!
      END SUBROUTINE CNDET_MCLR
