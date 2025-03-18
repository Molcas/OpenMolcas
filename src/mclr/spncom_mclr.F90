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
      SUBROUTINE SPNCOM_MCLR(iwork,NOPEN,MS2,NDET,IABDET,               &
     &                  IABUPP,IFLAG,PSSIGN,IPRCSF)
!
! Combinations of nopen unpaired electrons.Required
! spin projection MS2/2.
! JO 21-7-84
!    IFLAG = 1 : Only combinations ( in IABDET )
!    IFLAG = 2 : combinations and upper dets
!    IFLAG = 3 : Only upper dets
! A few revisions october 1988
! Upper dets added feb 1989
! Changed to combinations June 1992
!
! If PSSIGN differs from 0, spin combinations are assumed.
! we select as the unique determinants those with first electron
! having alpha spin
!
      Implicit None
      INTEGER IWORK(*)
      INTEGER NOPEN,MS2,NDET
      INTEGER IABDET(NOPEN,*),IABUPP(NOPEN,*)
      REAL*8 PSSIGN
      INTEGER IFLAG,IPRCSF

!     local variables
      INTEGER ADD
      INTEGER NUPPER,MX,I,J,NALPHA,MS2L,lUPPER,IEL
!
! LENGTH OF IWORK MUST BE AT LEAST NOPEN
!
      NDET=0
      NUPPER = 0
!
! Determinants are considered as binary numbers,1=alpha,0=beta
!
      MX=2 ** NOPEN
      IWORK(1:NOPEN+1) = 0
! Loop over all possible binary numbers
      DO 200 I=1,MX
!.. 1 : NEXT BINARY NUMBER
        ADD=1
        J=0
  190   CONTINUE
        J=J+1
        IF(IWORK(J).EQ.1) THEN
          IWORK(J)=0
        ELSE
          IWORK(J)=1
          ADD=0
        END IF
        IF( ADD .EQ. 1 ) GOTO 190
!
!.. 2 :  CORRECT SPIN PROJECTION ?
        NALPHA=0
        DO 180 J=1,NOPEN
          NALPHA=NALPHA+IWORK(J)
  180   CONTINUE
!
        IF(2*NALPHA-NOPEN.EQ.MS2.AND.                                   &
     &    .NOT.(PSSIGN.NE.0.0D0 .AND. IWORK(1).EQ.0)) THEN
          IF (IFLAG .LT. 3 ) THEN
            NDET=NDET+1
            IABDET(:,NDET) = IWORK(1:NOPEN)
          END IF
          IF (IFLAG .GT. 1 ) THEN
! UPPER DET ?
            MS2L = 0
            LUPPER = 1
!
            DO 10 IEL = 1,NOPEN
              IF (IWORK(IEL).EQ.1) THEN
                 MS2L = MS2L + 1
              ELSE
                 MS2L = MS2L - 1
              END IF
              IF( MS2L .LT. 0 ) LUPPER = 0
   10       CONTINUE
            IF( LUPPER .EQ. 1 ) THEN
              NUPPER = NUPPER + 1
              IABUPP(:,NUPPER) = IWORK(1:NOPEN)
            END IF
          END IF
        END  IF
!
  200 CONTINUE
!
!     XMSD2=DBLE(MS2)/2
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRCSF)
      END SUBROUTINE SPNCOM_MCLR
