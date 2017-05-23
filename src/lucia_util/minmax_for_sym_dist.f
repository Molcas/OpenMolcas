************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1997,1998, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDIST)
*
* A combination of NIGRP groups are given (IGRP)
*. Find MIN and MAX for symmetry in each group
*
* Jeppe Olsen, September 1997
*              April 1998     From  MINMAX_SM_GP
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Include blocks
#include "mxpdim.fh"
#include "strbas.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "csm.fh"
#include "WrkSpc.fh"
*. Input
      DIMENSION IGRP(NIGRP)
*.Output
      DIMENSION MNVAL(NIGRP),MXVAL(NIGRP)
*. Local scratch
C-jwk-cleanup      DIMENSION LSMGP(MXPOBS,MXPNGAS)
*
      NTEST = 0000
      IF(NTEST.GE.100) WRITE(6,*) ' >> Entering MINMAX_... <<'
*
      DO JGRP = 1, NIGRP
        MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
        MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
      END DO

*. Number of strings per sym and group
C     DO JGRP = 1, NIGRP
C       CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(JGRP)-1)*NSMST+1,
C    &               NSMST,LSMGP(1,JGRP))
C     END DO
C     IF(NTEST.GE.1000) THEN
C       WRITE(6,*) ' LSMGP '
C       CALL IWRTMA(LSMGP,NSMST,NIGRP,MXPOBS,NIGRP)
C     END IF
C. Max and min sym in each group
C     DO JGRP = 1, NIGRP
*
C       IMAX = 1
C       DO ISM=1, NSMST
C         IF(LSMGP(ISM,JGRP).GT.0) IMAX = ISM
C       END DO
C       MXVAL(JGRP) = IMAX
*
C       IMIN = NSMST
C       DO ISM = NSMST,1,-1
C        IF(LSMGP(ISM,JGRP).GT.0) IMIN = ISM
C       END DO
C       MNVAL(JGRP) = IMIN
C     END DO
*. Total number of symmetry distributions
      NDIST = 1
      DO JGRP = 1, NIGRP
        NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
      END DO
*
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
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving MINMAX_... <<'
*
      RETURN
      END
