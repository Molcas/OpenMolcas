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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE TRDNS2O(IVEC,JVEC,DPT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "sigma.fh"
      DIMENSION DPT2(*)
      DIMENSION IFCOUP(13,13)
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

C Remove this after debugging:
C     WRITE(*,*)' TRDNS2O Warning: Inactive-Active, Active-'//
C    &  'Secondary and Inactive-Secondary parts'
C     WRITE(*,*)' of CASPT2 density matrix contribution to '//
C    & '2nd order in perturbation'
C     WRITE(*,*)' theory are presently not properly debugged.'
C     RETURN
      CALL QENTER('TRDNS2O')

C Enter coupling cases for non-diagonal blocks:
      DO I=1,NCASES
       DO J=1,NCASES
        IFCOUP(I,J)=0
       END DO
      END DO
      IFCOUP(2,1)=1
      IFCOUP(3,1)=2
      IFCOUP(5,1)=3
      IFCOUP(6,1)=4
      IFCOUP(7,1)=5
      IFCOUP(6,2)=6
      IFCOUP(7,3)=7
      IFCOUP(5,4)=8
      IFCOUP(8,4)=9
      IFCOUP(9,4)=10
      IFCOUP(10,4)=11
      IFCOUP(11,4)=12
      IFCOUP(6,5)=13
      IFCOUP(7,5)=14
      IFCOUP(10,5)=15
      IFCOUP(11,5)=16
      IFCOUP(12,5)=23
      IFCOUP(13,5)=24
      IFCOUP(12,6)=17
      IFCOUP(13,7)=18
      IFCOUP(10,8)=19
      IFCOUP(11,9)=20
      IFCOUP(12,10)=21
      IFCOUP(13,11)=22

C If the G1 correction to the Fock matrix is used, then the
C inactive/virtual coupling elements (which are non-zero for the
C case of average CASSCF) cannot be used in the CASPT2 equations.
      IF(FOCKTYPE.EQ.'G1      ' .AND. (.NOT. G1SECIN)) THEN
        IFCOUP(12,5)=0
        IFCOUP(13,5)=0
      END IF

C Add to the block-off-diag parts of transition density matrix,
C    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>.
C i.e. inact/act, inact/virt and act/virt submatrices only,
C where IVEC, JVEC stands for the 1st-order perturbed CASPT2
C wave functions stored as vectors nr IVEC, JVEC on LUSOLV.

      NDPT2=0
      DO ISYM=1,NSYM
        NDPT2=NDPT2+NORB(ISYM)**2
      END DO

C Loop over ordering: First, <IVEC|...|JVEC>, then reverse.
C For each order, compute the upper-triangular blocks.
C For true density matrices, use symmetry of D-matrix.

C Transform to standard representation, contravariant form.
      CALL PTRTOC(0,IVEC,IVEC)
      IF(IVEC.NE.JVEC) CALL PTRTOC(0,JVEC,JVEC)
      NLOOP=2
      IF(IVEC.EQ.JVEC) NLOOP=1
      DO 1000 ILOOP=1,NLOOP
        IF(ILOOP.EQ.1) THEN
          IBRA=IVEC
          IKET=JVEC
        ELSE
          IBRA=JVEC
          IKET=IVEC
        END IF

C Loop over types and symmetry block of VEC1 vector:
      DO 400 ICASE1=1,13
        DO 400 ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) GOTO 400
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          NVEC1=NIS1*NAS1
          IF(NVEC1.EQ.0) GOTO 400
C Form VEC1 from the BRA vector, transformed to covariant form.
          CALL RHS_ALLO(NAS1,NIS1,LVEC1)
          CALL RHS_SCAL(NAS1,NIS1,LVEC1,0.0D0)
          IF(ICASE1.LE.11) THEN
           CALL RHS_ALLO(NAS1,NIS1,LSCR)
           CALL RHS_READ(NAS1,NIS1,LSCR,ICASE1,ISYM1,IBRA)
           CALL RHS_STRANS (NAS1,NIS1,1.0D0,LSCR,LVEC1,ICASE1,ISYM1)
           CALL RHS_FREE(NAS1,NIS1,LSCR)
          ELSE
           CALL RHS_READ(NAS1,NIS1,LVEC1,ICASE1,ISYM1,IBRA)
          END IF
C Form WEC1 from VEC1, if needed.
          NWEC1=0
          LWEC1=1
          FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
          IF(ICASE1.EQ.1) NWEC1=NASH(ISYM1)*NISH(ISYM1)
          IF(ICASE1.EQ.4) NWEC1=NASH(ISYM1)*NSSH(ISYM1)
          IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) NWEC1=NIS1
          IF(NWEC1.GT.0) THEN
            CALL GETMEM('WEC1','ALLO','REAL',LWEC1,NWEC1)
            CALL DCOPY_(NWEC1,[0.0D0],0,WORK(LWEC1),1)
            IMLTOP=1
#ifdef _MOLCAS_MPP_
            IF (IS_REAL_PAR()) THEN
                CALL GETMEM('TMP1','ALLO','REAL',LTMP1,NVEC1)
                CALL RHS_GET(NAS1,NIS1,LVEC1,WORK(LTMP1))
                IF(ICASE1.EQ.1) THEN
                  CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LTMP1),
     &                    WORK(LWEC1))
                ELSE IF(ICASE1.EQ.4) THEN
                  CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LTMP1),
     &                    WORK(LWEC1))
                ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
                  CALL SPEC1D(IMLTOP,FACT,WORK(LTMP1),WORK(LWEC1))
                END IF
                CALL GETMEM('TMP1','FREE','REAL',LTMP1,NVEC1)
            ELSE
              IF(ICASE1.EQ.1) THEN
                CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LVEC1),
     &                    WORK(LWEC1))
              ELSE IF(ICASE1.EQ.4) THEN
                CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LVEC1),
     &                    WORK(LWEC1))
              ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
                CALL SPEC1D(IMLTOP,FACT,WORK(LVEC1),WORK(LWEC1))
              END IF
            END IF
#else
            IF(ICASE1.EQ.1) THEN
              CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LVEC1),
     &                    WORK(LWEC1))
            ELSE IF(ICASE1.EQ.4) THEN
              CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LVEC1),
     &                    WORK(LWEC1))
            ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
              CALL SPEC1D(IMLTOP,FACT,WORK(LVEC1),WORK(LWEC1))
            END IF
#endif
          END IF
C Note: WEC1 is identical to <IBRA| E(p,q) |0> for the cases
C (p,q)=(t,i), (a,t), and (a,i), resp.
          DO 300 ICASE2=ICASE1+1,13
            IF(IFCOUP(ICASE2,ICASE1).EQ.0) GOTO 300
            DO 200 ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) GOTO 200
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NVEC2=NIS2*NAS2
              IF(NVEC2.EQ.0) GOTO 200
              CALL RHS_ALLO(NAS2,NIS2,LVEC2)
              CALL RHS_READ(NAS2,NIS2,LVEC2,ICASE2,ISYM2,IKET)
#ifdef _MOLCAS_MPP_
              IF (IS_REAL_PAR()) THEN
                  CALL GETMEM('TMP1','ALLO','REAL',LTMP1,NVEC1)
                  CALL GETMEM('TMP2','ALLO','REAL',LTMP2,NVEC2)
                  CALL RHS_GET(NAS1,NIS1,LVEC1,WORK(LTMP1))
                  CALL RHS_GET(NAS2,NIS2,LVEC2,WORK(LTMP2))
                  CALL OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,
     &                  WORK(LWEC1),WORK(LTMP1),DPT2,WORK(LTMP2),
     &                  iWORK(LLISTS))
                  CALL GETMEM('TMP1','FREE','REAL',LTMP1,NVEC1)
                  CALL GETMEM('TMP2','FREE','REAL',LTMP2,NVEC2)
              ELSE
                CALL OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,
     &                WORK(LWEC1),WORK(LVEC1),DPT2,WORK(LVEC2),
     &                iWORK(LLISTS))
              END IF
#else
              CALL OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,
     &              WORK(LWEC1),WORK(LVEC1),DPT2,WORK(LVEC2),
     &              iWORK(LLISTS))
#endif
              CALL RHS_FREE(NAS2,NIS2,LVEC2)
 200        CONTINUE
 300      CONTINUE
          CALL RHS_FREE(NAS1,NIS1,LVEC1)
          IF(NWEC1.GT.0)
     &         CALL GETMEM('WEC1','FREE','REAL',LWEC1,NWEC1)
 400  CONTINUE

      CALL GADSUM(DPT2,NDPT2)

      IF(IVEC.NE.JVEC) THEN
C Transpose the density matrix.
        CALL GETMEM('SCR','ALLO','REAL',LSCR,NDPT2)
        ISTA=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          CALL TRNSPS(NO,NO,DPT2(ISTA),WORK(LSCR))
          CALL DCOPY_(NO*NO,WORK(LSCR),1,DPT2(ISTA),1)
          ISTA=ISTA+NO**2
        END DO
        CALL GETMEM('SCR','FREE','REAL',LSCR,NDPT2)
      END IF
 1000 CONTINUE

C Transform vectors back to eigenbasis of H0(diag).
      CALL PTRTOSR(1,IVEC,IVEC)
      IF(IVEC.NE.JVEC) CALL PTRTOSR(1,JVEC,JVEC)
      IF(IVEC.EQ.JVEC) THEN
C Fill in lower-triangular block elements by symmetry.
        IDOFF=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NO=NORB(ISYM)
          DO IP=1,NI+NA
            DO IQ=NI+1,NO
              IDPQ=IDOFF+IP+NO*(IQ-1)
              IDQP=IDOFF+IQ+NO*(IP-1)
              DPT2(IDQP)=DPT2(IDPQ)
            END DO
          END DO
          IDOFF=IDOFF+NO**2
        END DO
      END IF

      CALL QEXIT('TRDNS2O')
      RETURN
      END
