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
      SUBROUTINE TRDNS2O(IVEC,JVEC,DPT2,NDPT2,SCAL)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: LISTS
      use EQSOLV, only: IfCoup
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: FockType, G1SecIn, nActEl, nSym, nOrb,
     &                         nInDep, nISup, nASup, nIsh, nSsh, nAsh

      IMPLICIT None

      Integer(kind=iwp), Intent(out):: NDPT2
      Integer(kind=iwp), Intent(in):: IVEC, JVEC
      real(kind=wp), intent(in):: SCAL
      real(kind=wp), intent(inout):: DPT2(*)

      real(kind=wp), Allocatable:: WEC1(:), SCR(:)
#ifdef _MOLCAS_MPP_
      real(kind=wp), Allocatable:: TMP1(:), TMP2(:)
#endif
      real(kind=wp) Fact
      Integer(kind=iwp) iCase1, iCase2, idoff, idpq, idqp, iLoop,
     &                  iMltOp, ip, iq, iSta, iSym, iSym1, iSym2,
     &                  lScr, lScr2, lVec1, lVec2, na, nas1, nas2,
     &                  ni, nis1, nis2, nLoop, no, nVec1, nVec2,
     &                  nWec1


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
      DO ILOOP=1,NLOOP
        ! IF(ILOOP.EQ.1) THEN
        !   IBRA=IVEC
        !   IKET=JVEC
        ! ELSE
        !   IBRA=JVEC
        !   IKET=IVEC
        ! END IF

C Loop over types and symmetry block of VEC1 vector:
      DO ICASE1=1,13
        DO ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) Cycle
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          NVEC1=NIS1*NAS1
          IF(NVEC1.EQ.0) Cycle
C Form VEC1 from the BRA vector, transformed to covariant form.
          CALL RHS_ALLO(NAS1,NIS1,LVEC1)
          CALL RHS_SCAL(NAS1,NIS1,LVEC1,Zero)
          IF(ICASE1.LE.11) THEN
           CALL RHS_ALLO(NAS1,NIS1,LSCR)
           CALL RHS_READ(NAS1,NIS1,LSCR,ICASE1,ISYM1,IVEC) !! IBRA)
           IF (IVEC.NE.JVEC.AND.ILOOP.EQ.1) THEN
            IF (SCAL.ne.One) CALL RHS_SCAL(NAS1,NIS1,LSCR,SCAL)
            CALL RHS_ALLO(NAS1,NIS1,LSCR2)
            CALL RHS_READ(NAS1,NIS1,LSCR2,ICASE1,ISYM1,JVEC)
            CALL RHS_DAXPY(NAS1,NIS1,One,LSCR2,LSCR)
            CALL RHS_FREE(LSCR2)
           END IF
           CALL RHS_STRANS (NAS1,NIS1,One,LSCR,LVEC1,ICASE1,ISYM1)
           CALL RHS_FREE(LSCR)
          ELSE
           CALL RHS_READ(NAS1,NIS1,LVEC1,ICASE1,ISYM1,IVEC) !! IBRA)
           IF (IVEC.NE.JVEC.AND.ILOOP.EQ.1) THEN
            IF (SCAL.ne.One) CALL RHS_SCAL(NAS1,NIS1,LVEC1,SCAL)
            CALL RHS_ALLO(NAS1,NIS1,LSCR2)
            CALL RHS_READ(NAS1,NIS1,LSCR2,ICASE1,ISYM1,JVEC)
            CALL RHS_DAXPY(NAS1,NIS1,One,LSCR2,LVEC1)
            CALL RHS_FREE(LSCR2)
           END IF
          END IF
C Form WEC1 from VEC1, if needed.
          NWEC1=0
          FACT=One/(DBLE(MAX(1,NACTEL)))
          IF(ICASE1.EQ.1) NWEC1=NASH(ISYM1)*NISH(ISYM1)
          IF(ICASE1.EQ.4) NWEC1=NASH(ISYM1)*NSSH(ISYM1)
          IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) NWEC1=NIS1
          IF(NWEC1.GT.0) THEN
            CALL mma_allocate(WEC1,NWEC1,Label='WEC1')
            WEC1(:)=Zero
            IMLTOP=1
#ifdef _MOLCAS_MPP_
            IF (IS_REAL_PAR()) THEN
                CALL mma_allocate(TMP1,NVEC1,Label='TMP1')
                CALL RHS_GET(NAS1,NIS1,LVEC1,TMP1)
                IF(ICASE1.EQ.1) THEN
                  CALL SPEC1A(IMLTOP,FACT,ISYM1,TMP1,WEC1)
                ELSE IF(ICASE1.EQ.4) THEN
                  CALL SPEC1C(IMLTOP,FACT,ISYM1,TMP1,WEC1)
                ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
                  CALL SPEC1D(IMLTOP,FACT,TMP1,WEC1)
                END IF
                CALL mma_deallocate(TMP1)
            ELSE
#endif
              IF(ICASE1.EQ.1) THEN
                CALL SPEC1A(IMLTOP,FACT,ISYM1,
     &                      GA_Arrays(LVEC1)%A,WEC1)
              ELSE IF(ICASE1.EQ.4) THEN
                CALL SPEC1C(IMLTOP,FACT,ISYM1,
     &                      GA_Arrays(LVEC1)%A,WEC1)
              ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
                CALL SPEC1D(IMLTOP,FACT,GA_Arrays(LVEC1)%A,WEC1)
              END IF
#ifdef _MOLCAS_MPP_
            END IF
#endif
          ELSE
            CALL mma_allocate(WEC1,1,Label='WEC1')
          END IF
C Note: WEC1 is identical to <IBRA| E(p,q) |0> for the cases
C (p,q)=(t,i), (a,t), and (a,i), resp.
          DO ICASE2=ICASE1+1,13
            IF(IFCOUP(ICASE2,ICASE1).EQ.0) Cycle
            DO ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) Cycle
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NVEC2=NIS2*NAS2
              IF(NVEC2.EQ.0) Cycle
              CALL RHS_ALLO(NAS2,NIS2,LVEC2)
              CALL RHS_READ(NAS2,NIS2,LVEC2,ICASE2,ISYM2,IVEC) !! IKET)
              IF (IVEC.NE.JVEC.AND.ILOOP.EQ.2) THEN
               IF (SCAL.ne.One) CALL RHS_SCAL(NAS2,NIS2,LVEC2,SCAL)
               CALL RHS_ALLO(NAS2,NIS2,LSCR2)
               CALL RHS_READ(NAS2,NIS2,LSCR2,ICASE2,ISYM2,JVEC)
               CALL RHS_DAXPY(NAS2,NIS2,One,LSCR2,LVEC2)
               CALL RHS_FREE(LSCR2)
              END IF
#ifdef _MOLCAS_MPP_
              IF (IS_REAL_PAR()) THEN
                  CALL mma_allocate(TMP1,NVEC1,Label='TMP1')
                  CALL mma_allocate(TMP2,NVEC2,Label='TMP2')
                  CALL RHS_GET(NAS1,NIS1,LVEC1,TMP1)
                  CALL RHS_GET(NAS2,NIS2,LVEC2,TMP2)
                  CALL OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,
     &                        WEC1,TMP1,DPT2,TMP2,LISTS)
                  CALL mma_deallocate(TMP1)
                  CALL mma_deallocate(TMP2)
              ELSE
#endif
                CALL OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,
     &                      WEC1,GA_Arrays(LVEC1)%A,DPT2,
     &                           GA_Arrays(LVEC2)%A,LISTS)
#ifdef _MOLCAS_MPP_
              END IF
#endif
              CALL RHS_FREE(LVEC2)
            End Do
          End Do
          CALL RHS_FREE(LVEC1)
          Call mma_deallocate(WEC1)
        End Do
      End Do

      End Do

      CALL GADSUM(DPT2,NDPT2)

      IF(IVEC.NE.JVEC) THEN
C Transpose the density matrix.
        CALL mma_allocate(SCR,NDPT2,Label='SCR')
        ISTA=1
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          CALL TRNSPS(NO,NO,DPT2(ISTA),SCR)
          CALL DCOPY_(NO*NO,SCR,1,DPT2(ISTA),1)
          ISTA=ISTA+NO**2
        END DO
        CALL mma_deallocate(SCR)
      END IF

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

      END SUBROUTINE TRDNS2O
