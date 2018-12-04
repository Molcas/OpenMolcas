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
* Copyright (C) 1997, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE TS_SYM_PNT2(   IGRP,  NIGRP, MAXVAL, MINVAL,   ISYM,
     &                          IPNT,   LPNT)
*
* Construct pointers to start of symmetrydistributions
* for supergroup of strings with given symmetry
*
* The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
* is given as
*     1
*     +  (ISM1-MINVAL(1))
*     +  (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     +  (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
*
* Where L is the last group of strings with nonvanishing occupation
*
* Jeppe Olsen, September 1997
*
* Version 2 : Uses IGRP and NIGRP to define supergroup
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
!      COMMON/BIGGY/WORK(MXPWRD)
#include "orbinp.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "strbas.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "csm.fh"
*. Specific Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
C-jwk-cleanup      INTEGER NELFGS(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
*. Output
      INTEGER MINVAL(*),MAXVAL(*),IPNT(*)
*
      NTEST = 00
*. Info on groups of strings in supergroup
      NGASL = 1
      DO IGAS = 1, NIGRP
       ITPFGS(IGAS) = IGRP(IGAS)
       IF(NELFGP(IGRP(IGAS)).GT.0) NGASL = IGAS
*. Number of strings per symmetry in each gasspace
C       CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
C    &               NNSTSGP(1,IGAS))
        CALL ICOPVE(NSTFSMGP(1,IGRP(IGAS)),NNSTSGP(1,IGAS),NSMST)
      END DO
*
C     NGASL = NIGRP
*
C     DO IGAS = 1, NIGRP
C       DO ISMST = 1, NSMST
C         IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
C       END DO
C       DO ISMST = NSMST,1,-1
C         IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
C       END DO
C     END DO
      DO IGAS = 1, NIGRP
        MINVAL(IGAS) = MINMAX_SM_GP(1,IGRP(IGAS))
        MAXVAL(IGAS) = MINMAX_SM_GP(2,IGRP(IGAS))
      END DO
*
      IF(NTEST.GE.1000) THEN
        write(6,*) 'NIGRP:', NIGRP
        WRITE(6,*)  ' MINVAL and MAXVAL '
        CALL IWRTMA(MINVAL,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MAXVAL,1,NIGRP,1,NIGRP)
      END IF

*. Total number of symmetry blocks that will be generated
      NBLKS = 1
      DO IGAS = 1, NGASL-1
       NBLKS = NBLKS*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
      END DO
      IF(NBLKS.GT.LPNT) THEN
        WRITE(6,*) ' Problem in TS_SYM_PNT'
        WRITE(6,*) ' Dimension of IPNT too small'
        WRITE(6,*) ' Actual and required length',NBLKS,LPNT
        WRITE(6,*)
        WRITE(6,*) ' I will Stop and wait for instructions'
*        STOP' TS_SYM_PNT too small'
        CALL SYSABENDMSG('lucia_util/ts_sym_pnt',
     &                    'Internal error',' ')
      END IF
*. Loop over symmetry blocks in standard order
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
C       ISTSMM1 = 1
C       DO IGAS = 1, NGASL -1
C         CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
C         ISTSMM1 = JSTSMM1
C       END DO
        ISTSMM1 = ISYMSTR(ISMFGS,NGASL-1)
*.  sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISYM)
        ISMFGS(NGASL) = ISMGSN
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL-1
          IOFF = IOFF + (ISMFGS(IGAS)-MINVAL(IGAS))*MULT
          MULT = MULT * (MAXVAL(IGAS)-MINVAL(IGAS)+1)
        END DO
*
        IPNT(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOFF, IPNT(IOFF) NSTRII ',
     &                 IOFF, IPNT(IOFF),NSTRII
        END IF
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TS_SYM_PNT'
        WRITE(6,*) ' Required total symmetry',ISYM
        WRITE(6,*) ' Number of symmetry blocks ', NBLKS
        WRITE(6,*)
        WRITE(6,*) ' Offset array  for symmetry blocks'
        CALL IWRTMA(IPNT,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END
