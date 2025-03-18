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
      SUBROUTINE DETSTR_MCLR(IDET,IASTR,IBSTR,NEL,NAEL,NBEL,NORB,       &
     &                       ISIGN,IWORK,IPRNT)

!
! A DETERMINANT,IDET,IS GIVEN AS A SET OF OCCUPIED SPIN ORBITALS,
! POSITIVE NUMBER INDICATES ALPHA ORBITAL AND NEGATIVE NUMBER
! INDICATES BETA ORBITAL .
!
! FIND CORRESPONDING ALPHA STRING AND BETA STRING ,
! AND DETERMINE SIGN NEEDED TO CHANGE DETERMINANT
! INTO PRODUCT OF ORDERED ALPHA STRING AND
! BETA STRING
!
! JEPPE OLSEN NOVEMBER 1988
!
      IMPLICIT NONE

      Integer NEL,NAEL,NBEL
      Integer IDET(NEL)
      Integer IASTR(NAEL),IBSTR(NBEL)
      Integer NORB,ISIGN
      Integer IWORK(*)
      Integer IPRNT
!
      INTEGER NTEST,IBEL,ITMP
!
!
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
!
! FIRST REORDER SPIN ORBITALS IN ASCENDING SEQUENCE
! THIS WILL AUTOMATICALLY SPLIT ALPHA AND BETASTRING
!
      CALL ORDSTR_MCLR(IDET,IWORK,NEL,ISIGN,IPRNT)
!
! ALPHA STRING IS LAST NAEL ORBITALS
      CALL iCOPY(NAEL,IWORK(NBEL+1),1,IASTR,1)
!
! BETA  STRING MUST BE COMPLETELY TURNED AROUND
      DO 10 IBEL = 1, NBEL
        IBSTR(IBEL) = -IWORK(NBEL+1-IBEL)
10    CONTINUE
! SIGN CHANGE FOR SWITCH OF BETA ORBITALS
      iTmp= NBEL*(NBEL+1)/2
      ISIGN = ISIGN * (-1) ** iTmp
!
      IF( NTEST.GE.200) THEN
        WRITE(6,*) ' INPUT DETERMINANT '
        CALL IWRTMA(IDET,1,NEL,1,NEL)
        WRITE(6,*) ' CORRESPONDING ALPHA STRING '
        CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
        WRITE(6,*) ' CORRESPONDING BETA STRING '
        CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
        WRITE(6,*) ' ISIGN FOR SWITCH ', ISIGN
      END IF

!      if(doDMRG.and.doMCLR)then ! yma
!        DO I=1,NEL
!          Write(117,1110,advance='no') IDET(I)
!        end do
!        Write(117,"(A,1X,I2)",advance='no')" SIGN",ISIGN
!1110  FORMAT(1X,I5)
!      end if

!
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NORB)
      END SUBROUTINE DETSTR_MCLR
