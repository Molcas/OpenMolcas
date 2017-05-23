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
* Copyright (C) 2015, Lasse Kragh Soerensen                            *
************************************************************************
      SUBROUTINE I_AM_SO_EXCITED(NBATCH,IBATCH,LBATCH,I1BATCH)
*
* Subroutine by Lasse from Oktober 2015
*
* Will give single excited states in from the desired GAS (or GAS's)
*
      IMPLICIT REAL*8(A-H,O-Z) ! I am so against this
*
#include "mxpdim.fh"
#include "gasstr.fh"
#include "cgas.fh"
*
* Input
      INTEGER IBATCH(8,*),LBATCH(*),I1BATCH(*)
* Scratch
      INTEGER IMAX_OCC(2,NGAS,2)
      INTEGER MAX_E_GAS_ALPHA(2,MXPSTT),MAX_E_GAS_BETA(2,MXPSTT)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Oh I am so excited '
        WRITE(6,*)
        WRITE(6,*) ' Number of GAS without maximum occupation = ',
     &               N_ELIMINATED_GAS
        WRITE(6,*)
        WRITE(6,*) ' GAS without maximum occupation '
        WRITE(6,*)
        DO I = 1, N_ELIMINATED_GAS
          WRITE(6,*) IELIMINATED_IN_GAS(I)
        END DO
      END IF
*
* First we need to find the GAS spaces for which we will eliminate
* the maximum occupation.
*
      IZERO = 0
      IMAX_OCC = IZERO
*
      DO JBATCH = 1,2 ! only alpha and beta
        DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
          ITYPE = IBATCH(JBATCH,IBLOCK)
          DO ISPGP = 1, NSPGPFTP(JBATCH)
            IOFF = IBSPGPFTP(JBATCH)
            DO IGAS = 1, NGAS
              IITYPE = ISPGPFTP(IGAS,IOFF-1+ISPGP)
              IEL = NELFGP(IITYPE)
              IF(IEL.GT.IMAX_OCC(JBATCH,IGAS,1))THEN
                IMAX_OCC(JBATCH,IGAS,1) = IEL
                IMAX_OCC(JBATCH,IGAS,2) = IITYPE
              END IF
            END DO
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        DO JBATCH = 1,2
          IF(JBATCH.EQ.1) THEN
            WRITE(6,*) ' Maximum number of alpha electrons in each GAS'
          ELSE
            WRITE(6,*) ' Maximum number of beta electrons in each GAS'
          END IF
          WRITE(6,*)
          WRITE(6,*) ' GAS, Electrons, Group '
          DO IGAS = 1,NGAS
            WRITE(6,*) IGAS,IMAX_OCC(JBATCH,IGAS,1),
     &                      IMAX_OCC(JBATCH,IGAS,2)
          END DO
          WRITE(6,*)
        END DO
      END IF
*
* Find which types contains the groups with a maximum number
* of alpha or beta electrons in a GAS
*
      NALPHA = 0
      NBETA = 0
      DO JBATCH = 1,2 ! only alpha and beta
        DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
          ITYPE = IBATCH(JBATCH,IBLOCK)
          DO ISPGP = 1, NSPGPFTP(JBATCH)
            IOFF = IBSPGPFTP(JBATCH)
            DO IGAS = 1, NGAS
              IITYPE = ISPGPFTP(IGAS,IOFF-1+ISPGP)
              IEL = NELFGP(IITYPE)
              IF(IEL.EQ.IMAX_OCC(JBATCH,IGAS,1))THEN
                IF(JBATCH.EQ.1) THEN
                  NALPHA = NALPHA + 1
                  MAX_E_GAS_ALPHA(1,NALPHA) = IGAS
                  MAX_E_GAS_ALPHA(2,NALPHA) = ISPGP
                ELSE
                  NBETA = NBETA + 1
                  MAX_E_GAS_BETA(1,NBETA) = IGAS
                  MAX_E_GAS_BETA(2,NBETA) = ISPGP
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Maximum number of alpha supergroups '//
     &             'that can be eliminated',NALPHA
        WRITE(6,*)
        WRITE(6,*) ' GAS Supergroup '
        WRITE(6,*)
        DO IGAS = 1, NGAS
          DO I = 1, NALPHA
            IF(MAX_E_GAS_ALPHA(1,I).EQ.IGAS) THEN
              WRITE(6,*) MAX_E_GAS_ALPHA(1,I),MAX_E_GAS_ALPHA(2,I)
            END IF
          END DO
        END DO
        WRITE(6,*)
        WRITE(6,*) 'Maximum number of beta supergroups '//
     &             'that can be eliminated',NBETA
        WRITE(6,*)
        WRITE(6,*) ' GAS Supergroup '
        WRITE(6,*)
        DO IGAS = 1, NGAS
          DO I = 1, NBETA
            IF(MAX_E_GAS_BETA(1,I).EQ.IGAS) THEN
              WRITE(6,*) MAX_E_GAS_BETA(1,I),MAX_E_GAS_BETA(2,I)
            END IF
          END DO
        END DO
        WRITE(6,*)
      END IF
*
* Now find the batches to possibly eliminate
*
      N_ELIMINATED_BATCHES = 0
*
      DO I = 1, N_ELIMINATED_GAS
        IGAS_ELIM = IELIMINATED_IN_GAS(I)
        DO JBATCH = 1, NBATCH
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            ITYPE_A = IBATCH(1,IBLOCK)
            ITYPE_B = IBATCH(2,IBLOCK)
* Will first check if it matches a beta type
            IMATCH_BETA = 0
            DO J = 1, NBETA
              IF(ITYPE_B.EQ.MAX_E_GAS_BETA(2,J).AND.
     &          IGAS_ELIM.EQ.MAX_E_GAS_BETA(1,J)) THEN
                IMATCH_BETA = 1
              ELSE
                CYCLE
              END IF
              IF(IMATCH_BETA.EQ.1) EXIT
            END DO
* Now check it also matches an alpha type (if the beta type is matched)
            IMATCH_ALPHA = 0
            IF(IMATCH_BETA.EQ.1) THEN
              DO J = 1, NALPHA
                IF(ITYPE_A.EQ.MAX_E_GAS_ALPHA(2,J).AND.
     &            IGAS_ELIM.EQ.MAX_E_GAS_ALPHA(1,J)) THEN
                  IMATCH_ALPHA = 1
                ELSE
                  CYCLE
                END IF
                IF(IMATCH_ALPHA.EQ.1) EXIT
              END DO
            END IF
            IF(IMATCH_BETA.EQ.1.AND.IMATCH_ALPHA.EQ.1) THEN
              N_ELIMINATED_BATCHES = N_ELIMINATED_BATCHES + 1
              I_AM_OUT(N_ELIMINATED_BATCHES) = JBATCH
            END IF
          END DO
        END DO
      END DO
*
      IF(N_ELIMINATED_BATCHES.GT.MXPSTT) THEN
        WRITE(6,*) ' Increase MXPSTT to ',N_ELIMINATED_BATCHES
        CALL SYSABENDMSG('lucia_util/i_am_so_excited',
     &                   'Dimension of I_AM_OUT is too small',
     &                   'Increase MXPSTT')
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of eliminated batches ',N_ELIMINATED_BATCHES
        WRITE(6,*)
        WRITE(6,*) ' The batches eliminated '
        WRITE(6,*)
        DO I = 1, N_ELIMINATED_BATCHES
          WRITE(6,*) I_AM_OUT(I)
        END DO
        WRITE(6,*)
      END IF
*
      END SUBROUTINE I_AM_SO_EXCITED
