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
! Copyright (C) 1987, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!***********************************************************************
      SUBROUTINE CNDET(ICONF,IPDET,NDET,NEL,NORB,NOP,NCL,               &
     &                 IDET,IPRINT)
!
!     AUTHOR:        J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
!     MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                    M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
!     PURPOSE:
!
!     A CONFIGURATION ICONF IN COMPRESSED FORM AND A SET OF
!     PROTOTYPE DETERMINANTS ,IPDET, ARE GIVEN .
!     CONSTRUCT THE CORRESPONDING DETERMINANTS IN CONTRACTED  FORM .
!
!     SUBROUTINE CALLS:
!
!     IWRTMA
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION ICONF(*   )
      DIMENSION IPDET(*       )
      DIMENSION IDET(NEL,*   )
!
! POSITIVE NUMBER  : ALPHA ORBITAL
! NEGATIVE NUMBER  : BETA  ORBITAL
!
!      IPRINT=40
!      write(*,*)"iprint eq 40 in CNDET, nclosed, nopen",iprint, ncl,nop
!      write(*,*)"iconf, nclosed, nopen",iconf(1:1+ncl-1)
!      if(nop.ne.0)then
!        write(*,*)"opened orbitals",iconf(1+ncl:ncl+nop)
!      end if

      IF( IPRINT.EQ.40 ) THEN
        IF(NCL .NE. 0 ) THEN
          Write(6,*) ' DOUBLE OCCUPIED ORBITALS '
          CALL IWRTMA(ICONF,1,NCL,1,NCL)
        END IF
        IF(NOP .NE. 0 ) THEN
          Write(6,*) ' OPEN ORBITALS '
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
!
          IF( IPDET(IADR    ) .EQ. 0 ) IDET(IADD+IOP,JDET) =            &
     &    - ICONF(NCL +IOP)
!
190     CONTINUE
200   CONTINUE
!
      IF( IPRINT.EQ.40 ) THEN
       Write(6,*) ' CONFIGURATION FROM DETCON '
       CALL IWRTMA(ICONF,1,NORB,1,NORB)
       Write(6,* ) ' PROTO TYPE DETERMINANTS '
       IF(NOP*NDET .GT. 0)                                              &
     & CALL IWRTMA(IPDET,NOP,NDET,NOP,NDET)
!
       IF(NEL*NDET .GT. 0 )                                             &
     & Write(6,*) ' CORRESPONDING DETERMINANTS '
       CALL IWRTMA(IDET,NEL,NDET,NEL,NDET)
!
      END IF

      IPRINT=0
!
      RETURN
      END
