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
* Copyright (C) 2001, Jeppe Olsen                                      *
*               2001, Giovanni Li Manni                                *
*               2001, Dongxia Ma                                       *
************************************************************************
      SUBROUTINE NSTPTP_GAS(NGAS,ISPGRP,NSTSGP,NSMST,
     &                      NSTSSPGP,IGRP,MXNSTR,
     &                      NSMCLS,NSMCLSE,NSMCLSE1)
*
* Find number of strings per symmetry for the supergroup defined
* by the groups of ISPGRP. The obtained number of strings per sym
* is stored in NSTSSPGP(*,IGRP)
*
* Jeppe Olsen, Giovanni Li Manni, Dongxia Ma
* Dicember 2011 - the old version is too slow for many GAS spaces
*                 (new version simpler and quicker)
*
*. Also delivered:
*
* NSMCLS : MAX Number of symmetry classes for given supergroup,
*          i.e. number of combinations of symmetries of groups
*          containing strings
* NSMCLSE : Number of symmetry classes for given supergroup
*          obtained by restricting allowed symmetries in
*          a given group by a max and min.
* NSMCLSE1 : As NSMCLSE, but the symmetry of the last active
*            orbital space where there is more than one symmetry
*            is left out
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ISPGRP(NGAS),NSTSGP(NSMST,*)
*. Input and Output (column IGRP updated)
      DIMENSION NSTSSPGP(NSMST,IGRP)
#include "mxpdim.fh"
#include "multd2h.fh"
*. Scratch
      INTEGER ISM,MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER MSM1(MXPNSMST),MSM2(MXPNSMST)
      INTEGER ISM1(MXPNSMST),ISM2(MXPNSMST)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NSTPTP_GAS is speaking '
        WRITE(6,*) ' ======================='
*
        WRITE(6,*) ' Supergroup in action '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
      END IF
*
*. The NSMCLS* parameters
*
*. Max and min for each GASpace


      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
      END DO
*. Last space with more than one symmetry
      NGASL = 1
      DO IGAS = 1, NGAS
        IF(MXSM(IGAS).NE.MNSM(IGAS)) NGASL = IGAS
      END DO
*. NSMCLSE
      NSMCLSE = 1
      DO IGAS = 1, NGAS
        NSMCLSE = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE
      END DO
*. NSMCLSE1
      NSMCLSE1 = 1
      DO IGAS = 1, NGASL-1
        NSMCLSE1 = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE1
      END DO
      IZERO = 0
      DO IGAS = 1, NGAS
*. In ISM1, the number of strings per symmetry for the first
*  IGAS-1 spaces are given, obtain in ISM2 the number of strings per sym
*  for the first IGAS spaces
*. Also: in MSM1, MSM2, counts the number of nontrivial combinations per
*  sym
        IF(IGAS.EQ.1) THEN
*. ISM1: The number of strings per symmetry for zero electrons
         CALL ISETVC(ISM1,IZERO,NSMST)
         ISM1(1) = 1
         CALL ISETVC(MSM1,IZERO,NSMST)
         MSM1(1) = 1
        ELSE
*. copy from the ISM2 obtained for preceeding IGAS
         CALL ICOPVE(ISM2,ISM1,NSMST)
         CALL ICOPVE(MSM2,MSM1,NSMST)
        END IF
        CALL ISETVC(ISM2,IZERO,NSMST)
        CALL ISETVC(MSM2,IZERO,NSMST)
        DO ISM_IGASM1 = 1, NSMST
         DO ISM_IGAS = 1, NSMST
           ISM = MULTD2H(ISM_IGASM1,ISM_IGAS)
           ISM2(ISM) = ISM2(ISM) +
     &     ISM1(ISM_IGASM1)*NSTSGP(ISM_IGAS,ISPGRP(IGAS))
           IF(ISM1(ISM_IGASM1)*NSTSGP(ISM_IGAS,ISPGRP(IGAS)).NE.0)
     &     MSM2(ISM) = MSM2(ISM) + MSM1(ISM_IGASM1)
         END DO
        END DO
      END DO !loop over IGAS
      CALL ICOPVE(ISM2,NSTSSPGP(1,IGRP),NSMST)
*
      MXNSTR = 0
      NSMCLS = 0
      DO ISTRSM = 1, NSMST
        MXNSTR = MAX(MXNSTR,NSTSSPGP(ISTRSM,IGRP))
        NSMCLS = MAX(NSMCLS,MSM2(ISTRSM))
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     &  ' Number of strings per symmetry for supergroup',IGRP
        CALL IWRTMA10(NSTSSPGP(1,IGRP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Largest number of strings of given sym ',MXNSTR
*
        WRITE(6,'(A,3(2X,I8))') ' NSMCLS,NSMCLSE,NSMCLSE1=',
     &                       NSMCLS,NSMCLSE,NSMCLSE1
      END IF
      RETURN
      END
