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
! Copyright (C) 1997,1998, Jeppe Olsen                                 *
!***********************************************************************
      SUBROUTINE MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDIST)
!
! A combination of NIGRP groups are given (IGRP)
!. Find MIN and MAX for symmetry in each group
!
! Jeppe Olsen, September 1997
!              April 1998     From  MINMAX_SM_GP
!
!
      use lucia_data, only: MINMAX_SM_GP
      IMPLICIT NONE
!. Input
      INTEGER NIGRP,NDIST
      INTEGER IGRP(NIGRP)
!.Output
      INTEGER MNVAL(NIGRP),MXVAL(NIGRP)
!. Local scratch
      INTEGER NTEST,JGRP
!
      NTEST = 0000
      IF(NTEST.GE.100) WRITE(6,*) ' >> Entering MINMAX_... <<'
!
      DO JGRP = 1, NIGRP
        MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
        MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
      END DO

!. Number of strings per sym and group
!     DO JGRP = 1, NIGRP
!       CALL ICOPVE2(NSTSGP(1)%I,(IGRP(JGRP)-1)*NSMST+1,
!    &               NSMST,LSMGP(1,JGRP))
!     END DO
!     IF(NTEST.GE.1000) THEN
!       WRITE(6,*) ' LSMGP '
!       CALL IWRTMA(LSMGP,NSMST,NIGRP,MXPOBS,NIGRP)
!     END IF
!. Max and min sym in each group
!     DO JGRP = 1, NIGRP
!
!       IMAX = 1
!       DO ISM=1, NSMST
!         IF(LSMGP(ISM,JGRP).GT.0) IMAX = ISM
!       END DO
!       MXVAL(JGRP) = IMAX
!
!       IMIN = NSMST
!       DO ISM = NSMST,1,-1
!        IF(LSMGP(ISM,JGRP).GT.0) IMIN = ISM
!       END DO
!       MNVAL(JGRP) = IMIN
!     END DO
!. Total number of symmetry distributions
      NDIST = 1
      DO JGRP = 1, NIGRP
        NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
      END DO
!
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Group combination : '
        WRITE(6,'(5X,10I3)') (IGRP(JGRP),JGRP=1, NIGRP)
        WRITE(6,*)
        WRITE(6,*) ' Group Minsym Maxsym'
        WRITE(6,*) ' ==================='
        DO JGRP = 1, NIGRP
          WRITE(6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Total number of distributions', NDIST
      END IF
!
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving MINMAX_... <<'
!
      END SUBROUTINE MINMAX_FOR_SYM_DIST
