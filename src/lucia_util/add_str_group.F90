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
! Copyright (C) Jeppe Olsen                                            *
!***********************************************************************
      SUBROUTINE ADD_STR_GROUP( NSTADD, IOFADD, ISTADD,   NSTB,   NSTA, &
     &                         ISTRING,  IELOF, NELADD, NELTOT)
!
! Part of assembling strings in individual types to
! super group of strings
!
!. Copying strings belonging to a given type to supergroup of strings
!
! Jeppe Olsen, for once improving performance of LUCIA
!
!.Input
! =====
! NSTADD : Number of strings to be added
! IOFADD : First string to be added
! ISTADD : Strings to be added
! NSTB   : Number of strings belonging to lower gasspaces
! NSTA   : Number of strings belonging to higher gasspaces
! ISTRING: Supergroup of strings under construction
! IELOF  : Place of first electron to be added
! NELADD : Number of electrons to be added
! NELTOT : Total number of electrons
!
      IMPLICIT REAL*8(A-H,O-Z)
!. Input
!     DIMENSION ISTADD(NELADD,*)
      DIMENSION ISTADD(*)
!. Input and output
!     DIMENSION ISTRING(NELTOT,*)
      DIMENSION ISTRING(*)
!
!?    WRITE(6,*) '  Inside ADD ... '
      IF(NSTA.GT.1) THEN
        DO IISTR = 1,NSTADD
!. Address of A(1,IISTR,1)
!. A(I(after),Igas,I(before))
          IOFFY = (IOFADD-2+IISTR)*NELADD
          IOFF1 = (IISTR-1)*NSTA + 1
          IADD2 = NSTADD*NSTA
          IOFF2 = IOFF1 - IADD2
          DO ISTB = 1, NSTB
!. Address of A(1,IISTR,ISTB)
!           IOFF2 = IOFF1 + (ISTB-1)*NSTADD*NSTA
            IOFF2 = IOFF2 + IADD2
            IOFFX = IELOF-1+(IOFF2-2)*NELTOT
            DO ISTA = 1, NSTA
              IOFFX = IOFFX + NELTOT
              DO JEL = 1, NELADD
                ISTRING(JEL+IOFFX)                                      &
     &        = ISTADD(JEL+IOFFY)
!               ISTRING(IELOF-1+JEL,IOFF2-1+ISTA)
!    &        = ISTADD(JEL,IOFADD-1+IISTR)
              END DO
            END DO
          END DO
        END DO
      ELSE IF (NSTA .EQ. 1 ) THEN
!. Address of A(1,IISTR,1)
!. A(I(after),Igas,I(before))
        DO ISTB = 1, NSTB
          IOFF0 = (ISTB-1)*NSTADD
          IOFFY  = (IOFADD-2)*NELADD
          IOFFX = IELOF-1+(IOFF0-1)*NELTOT
          DO IISTR = 1,NSTADD
!. Address of A(1,IISTR,ISTB)
!           IOFF2 = IISTR  + IOFF0
!           IOFFX  = IELOF-1+(IOFF2-1)*NELTOT
            IOFFX  = IOFFX + NELTOT
            IOFFY  = IOFFY + NELADD
            DO JEL = 1, NELADD
              ISTRING(JEL+IOFFX)                                        &
     &      = ISTADD( JEL+IOFFY)
!             ISTRING(IELOF-1+JEL+(IOFF2-1)*NELTOT)
!    &      = ISTADD(JEL+(IOFADD-1+IISTR-1)*NELADD)
!             ISTRING(IELOF-1+JEL,IOFF2)
!    &      = ISTADD(JEL,IOFADD-1+IISTR)
            END DO
!?          WRITE(6,*) ' New string from ADD '
!?          CALL IWRTMA(ISTRING(IOFFX+1),1,NELADD,1,NELADD)
!?          WRITE(6,*) ' IOFFX, IOFFY, NELADD',
!?   &                   IOFFX, IOFFY, NELADD
          END DO
        END DO
      END IF
!
      RETURN
      END
