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
* Copyright (C) 1994,1999, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
* 1999: GEMINAL R12 ENABLED                  *
*--------------------------------------------*
      SUBROUTINE SIGMA_CASPT2(ALPHA,BETA,IVEC,JVEC)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use Fockof, only: FIT, FAI_Full, FIA_Full, FIT_Full, FTA_Full,
     &                  IOFFIT, IOFFIA, IOFFTA, FTI, FIA, FTI_Full,
     &                  FAI, FTA, FAT, FAT_Full
      use caspt2_global, only: FIFA, LISTS
      use stdalloc, only: mma_allocate, mma_deallocate
      use EQSOLV, only: IFCoup
      use Sigma_data, only: IFTEST, NFDXP, NFMV, NFR1, NFSCA
      use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array,
     &                   GA_Arrays
      use caspt2_module, only: CPUSGM, TIOSGM, FockType, G1SecIn, MaxIt,
     &                         nActEl, nCases, nSym, ThrShn, ThrShS,
     &                         nIsh, nAsh, nSsh, nOrb, nInDep, nISup,
     &                         nASup
      IMPLICIT None
      real(kind=wp), intent(in) :: ALPHA, BETA
      integer(kind=iwp), intent(in) :: IVEC, JVEC

      real(kind=wp), ALLOCATABLE:: SGM1(:), SGM2(:), D1(:), D2(:)
      real(kind=wp) CPU, CPU0, CPU1
      real(kind=wp) TIO, TIO0, TIO1
      real(kind=wp) Fact, XTST
      real(kind=wp), external:: RHS_DDot, DDot_
      integer(kind=iwp) iCASE1, iCase2, IfC, IFIFA, IMLTOP, ISYM,
     &                  ISYM1, ISYM2, lCX, lg_CX, lg_D2, lg_Sgm2,
     &                  lg_SgmX, lSgm2_Sta, lSgmX, lSgmX_Sta,
     &                  Max_MESG_Size, NA, NAS1, NAS2, nCX, ND1,
     &                  ND2, NFIA, NFIT, NFTA, NI, NIS1, NIS2, NO, NS,
     &                  nSgm1, nSgm2, nSgm2_Blk, nSgmX, nSgmX_blk

C Compute |JVEC> := BETA* |JVEC> + ALPHA* (H0-E0)* |IVEC>
C where the vectors are represented in transformed basis and
C are  stored at positions IVEC and JVEC on the LUSOLV unit.


#ifdef _DEBUGPRINT_
      WRITE(6,*)' Entering SIGMA.'
      WRITE(6,*)
     &' Compute |JVEC> := Beta*|JVEC> + Alpha*(H0-E0)|IVEC>'
      WRITE(6,'(1x,a,2f15.6)')'Alpha,Beta:',Alpha,Beta
      WRITE(6,'(1x,a,2i5)')'IVEC,JVEC:',IVEC,JVEC
#endif


C If the G1 correction to the Fock matrix is used, then the
C inactive/virtual coupling elements (which are non-zero for the
C case of average CASSCF) cannot be used in the CASPT2 equations.
      IF(FOCKTYPE.EQ.'G1      ' .AND. (.NOT. G1SECIN)) THEN
        IFCOUP(12,5)=0
        IFCOUP(13,5)=0
      END IF

      IFTEST=0
C Flop counts:
      NFSCA=0
      NFDXP=0
      NFMV =0
      NFR1 =0
C First compute diagonal block contributions:
CTEST      WRITE(6,*)' First, do it for (H0(diag)-E0).'
      CALL PSGMDIA(ALPHA,BETA,IVEC,JVEC)
      IF(ALPHA.EQ.Zero) Return
      IF(MAXIT.EQ.0) Return
CTEST      WRITE(6,*)
CTEST     & ' From now on, scaling with BETA is already done.'
CTEST      WRITE(6,*)' Test print  after SGMDIA call in SIGMA:'
CTEST      WRITE(6,*)' Should be zero, in first iteration.'
CTEST      CALL OVLPRT(JVEC,JVEC,OVLAPS)
C From now on, scaling with BETA is already done.

C Transform to standard representation:
      CALL PTRTOC(0,IVEC,IVEC)
      CALL PTRTOC(1,JVEC,JVEC)

C Set up non-diagonal blocks of Fock matrix:
C SVC: add transposed fock matrix blocks
      NFIT=0
      NFIA=0
      NFTA=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NS=NSSH(ISYM)
        IOFFIT(ISYM)=NFIT
        IOFFIA(ISYM)=NFIA
        IOFFTA(ISYM)=NFTA
        NFIT=NFIT+NA*NI
        NFIA=NFIA+NS*NI
        NFTA=NFTA+NS*NA
      END DO
      NFIT=NFIT+1 !?
      NFIA=NFIA+1 !?
      NFTA=NFTA+1 !?

      Call mma_allocate(FIT_Full,NFIT,Label='FIT_Full')
      Call mma_allocate(FTI_Full,NFIT,Label='FTI_Full')

      Call mma_allocate(FIA_Full,NFIA,Label='FIA_Full')
      Call mma_allocate(FAI_Full,NFIA,Label='FAI_Full')

      Call mma_allocate(FTA_Full,NFTA,Label='FTA_Full')
      Call mma_allocate(FAT_Full,NFTA,Label='FAT_Full')

      IFIFA=1
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NS=NSSH(ISYM)
        NO=NORB(ISYM)

        IF (NO > 0) THEN
        FIT(ISYM)%A(1:NA*NI) =>
     &     FIT_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)
        FTI(ISYM)%A(1:NA*NI) =>
     &     FTI_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)

        FIA(ISYM)%A(1:NS*NI) =>
     &     FIA_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)
        FAI(ISYM)%A(1:NS*NI) =>
     &     FAI_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)

        FTA(ISYM)%A(1:NS*NA) =>
     &     FTA_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)
        FAT(ISYM)%A(1:NS*NA) =>
     &     FAT_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)

        CALL FBLOCK(FIFA(IFIFA),NO,NI,NA,NS,
     &              FIT(ISYM)%A(:),FTI(ISYM)%A(:),
     &              FIA(ISYM)%A(:),FAI(ISYM)%A(:),
     &              FTA(ISYM)%A(:),FAT(ISYM)%A(:))

        IFIFA=IFIFA+(NO*(NO+1))/2
        END IF

      END DO

      CALL TIMING(CPU0,CPU,TIO0,TIO)
C Loop over types and symmetry block of sigma vector:
      DO ICASE1=1,11
        DO ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) Cycle
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          NSGM2=NIS1*NAS1
          IF(NSGM2.EQ.0) Cycle

          CALL mma_allocate(SGM2,NSGM2,Label='SGM2')
          SGM2(:)=Zero

          IF(ICASE1.EQ.1) THEN
            NSGM1=NASH(ISYM1)*NISH(ISYM1)
          ELSE IF(ICASE1.EQ.4) THEN
            NSGM1=NASH(ISYM1)*NSSH(ISYM1)
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            NSGM1=NIS1
          ELSE
            NSGM1=0
          END IF
          CALL mma_allocate(SGM1,MAX(1,NSGM1),LABEL='SGM1')
          SGM1(:)=Zero

          IMLTOP=0
          DO ICASE2=ICASE1+1,NCASES
            IFC=IFCOUP(ICASE2,ICASE1)
            IF(IFC.EQ.0) Cycle
            DO ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) Cycle
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NCX=NIS2*NAS2
              IF(NCX.EQ.0) Cycle

              CALL RHS_ALLO(NAS2,NIS2,lg_CX)
              CALL RHS_READ(NAS2,NIS2,lg_CX,ICASE2,ISYM2,IVEC)
C SVC: for case H (12,13) we can now pass the distributed array ID to
C the SGM subroutines
              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                LCX=lg_CX
                XTST=RHS_DDOT(NAS2,NIS2,lg_CX,lg_CX)
              ELSE
                LCX=Allocate_GA_Array(NCX,'CX')
                CALL RHS_GET(NAS2,NIS2,lg_CX,GA_Arrays(LCX)%A)
                CALL RHS_FREE(lg_CX)
                XTST=DDOT_(NCX,GA_Arrays(LCX)%A,1,
     &                         GA_Arrays(LCX)%A,1)
              END IF

              IF(XTST.GT.1.0D12) THEN
                WRITE(6,'(1x,a,6i10)')' SIGMA A. ICASE2,ISYM2:',
     &                                           ICASE2,ISYM2
                Call Crash()
              END IF

#ifdef _DEBUGPRINT_
              WRITE(6,*)' ISYM1,ICASE1:',ISYM1,ICASE1
              WRITE(6,*)' ISYM2,ICASE2:',ISYM2,ICASE2
              WRITE(6,*)' SIGMA calling SGM with IMLTOP=',IMLTOP
#endif
C Compute contribution SGM2 <- CX, and SGM1 <- CX  if any
              CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &                 SGM1,SGM2,LCX,LISTS)

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                CALL RHS_FREE(lg_CX)
              ELSE
                Call Deallocate_GA_Array(LCX)
              END IF

C Check for colossal values of SGM2 and SGM1
              XTST=DDOT_(NSGM2,SGM2,1,SGM2,1)
              IF(XTST.GT.1.0D12) THEN
                WRITE(6,'(1x,a,6i10)')' SIGMA B. ICASE1,ISYM1:',
     &                                           ICASE1,ISYM1
                WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
     &                                           ICASE2,ISYM2
                CALL Crash()
              END IF

              IF(NSGM1.GT.0) THEN
                XTST=DDOT_(NSGM1,SGM1,1,SGM1,1)
                IF(XTST.GT.1.0D12) THEN
                  WRITE(6,'(1x,a,6i10)')' SIGMA B2. ICASE1,ISYM1:',
     &                                              ICASE1,ISYM1
                  WRITE(6,'(1x,a,6i10)')'           ICASE2,ISYM2:',
     &                                              ICASE2,ISYM2
                  Call Crash()
                END IF
              END IF

            End Do
          End Do

C-SVC: sum the replicate arrays:
          MAX_MESG_SIZE = 2**27
          DO LSGM2_STA=1,NSGM2,MAX_MESG_SIZE
            NSGM2_BLK=MIN(MAX_MESG_SIZE,NSGM2-LSGM2_STA+1)
            CALL GADSUM(SGM2(LSGM2_STA:),NSGM2_BLK)
          END DO

          IF (NSGM1.GT.0) THEN
            CALL GADSUM(SGM1,NSGM1)
          END IF

C       XTST2=DDOT_(NSGM2,SGM2,1,SGM2,1)
C       XTST1=Zero
C       IF(NSGM1.GT.0)XTST1=DDOT_(NSGM1,SGM1,1,SGM1,1)
C       WRITE(6,'(1x,a,a,i2,2f16.6)')
C    & 'Contr. SGM2, SGM1, ',cases(icase1),isym1,xtst2,xtst1

C If there are 1-electron contributions, add them into the 2-el
C part (This requires a non-empty active space.)
          IF(NSGM1.GT.0) THEN
            FACT=One/(DBLE(MAX(1,NACTEL)))
            IF (ICASE1.EQ.1) THEN
              CALL SPEC1A(IMLTOP,FACT,ISYM1,SGM2,SGM1)
            ELSE IF(ICASE1.EQ.4) THEN
              CALL SPEC1C(IMLTOP,FACT,ISYM1,SGM2,SGM1)
            ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
              CALL SPEC1D(IMLTOP,FACT,SGM2,SGM1)
            END IF

            XTST=DDOT_(NSGM2,SGM2,1,SGM2,1)
            IF(XTST.GT.1.0D12) THEN
              WRITE(6,'(1x,a,6i10)')' SIGMA C. ICASE1,ISYM1:',
     &                                         ICASE1,ISYM1
              Call Crash()
            END IF

          END IF
          CALL mma_deallocate(SGM1)

C-SVC: no need for the replicate arrays any more, fall back to one array
          CALL RHS_ALLO (NAS1,NIS1,lg_SGM2)
          CALL RHS_PUT (NAS1,NIS1,lg_SGM2,SGM2)
          CALL mma_deallocate(SGM2)

C Add to sigma array. Multiply by S to  lower index.
          NSGMX=NSGM2
          CALL RHS_ALLO(NAS1,NIS1,lg_SGMX)
          CALL RHS_READ(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)

          XTST=RHS_DDOT(NAS1,NIS1,lg_SGMX,lg_SGMX)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA D. ICASE1,ISYM1:',ICASE1,ISYM1
            WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',ICASE2,ISYM2
            Call Crash()
          END IF

*         IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
            CALL RHS_STRANS(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX,
     &                      ICASE1,ISYM1)
*         ELSE
*           CALL RHS_DAXPY(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX)
*         END IF
          CALL RHS_FREE (lg_SGM2)

          XTST=RHS_DDOT(NAS1,NIS1,lg_SGMX,lg_SGMX)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA E. ICASE1,ISYM1:',ICASE1,ISYM1
            WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',ICASE2,ISYM2
            Call Crash()
          END IF

C Write SGMX to disk.
          CALL RHS_SAVE (NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
          CALL RHS_FREE (lg_SGMX)
        End Do
      End Do

      IMLTOP=1
C Loop over types and symmetry block of CX vector:
      DO ICASE1=1,11
        DO ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) Cycle
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          ND2=NIS1*NAS1
          IF(ND2.EQ.0) Cycle

          CALL RHS_ALLO (NAS1,NIS1,lg_D2)
          CALL RHS_SCAL (NAS1,NIS1,lg_D2,Zero)
C Contract S*CX to form D2. Also form D1 from D2, if needed.

          NCX=ND2
          CALL RHS_ALLO (NAS1,NIS1,lg_CX)
          CALL RHS_READ (NAS1,NIS1,lg_CX,ICASE1,ISYM1,IVEC)

          XTST=RHS_DDOT(NAS1,NIS1,lg_CX,lg_CX)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA F. ICASE1,ISYM1:',ICASE1,ISYM1
            Call Crash()
          END IF

          IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
           CALL RHS_STRANS (NAS1,NIS1,ALPHA,lg_CX,lg_D2,ICASE1,ISYM1)
          ELSE
           CALL RHS_DAXPY(NAS1,NIS1,ALPHA,lg_CX,lg_D2)
          END IF
          CALL RHS_FREE (lg_CX)

CPAM Sanity check:
          XTST=RHS_DDOT(NAS1,NIS1,lg_D2,lg_D2)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA G1 ICASE1,ISYM1:',ICASE1,ISYM1
            WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',ICASE2,ISYM2
            Call Crash()
          END IF

          CALL mma_allocate(D2,ND2,Label='D2')
          CALL RHS_GET (NAS1,NIS1,lg_D2,D2)
          CALL RHS_FREE (lg_D2)

          ND1=0
          IMLTOP=1
          FACT=One/(DBLE(MAX(1,NACTEL)))
          IF(ICASE1.EQ.1) THEN
            ND1=NASH(ISYM1)*NISH(ISYM1)
            IF(ND1.GT.0) THEN
              CALL mma_allocate(D1,ND1,Label='D1')
              D1(:)=Zero
              CALL SPEC1A(IMLTOP,FACT,ISYM1,D2,D1)
            END IF
          ELSE IF(ICASE1.EQ.4) THEN
            ND1=NASH(ISYM1)*NSSH(ISYM1)
            IF(ND1.GT.0) THEN
              CALL mma_allocate(D1,ND1,Label='D1')
              D1(:)=Zero
              CALL SPEC1C(IMLTOP,FACT,ISYM1,D2,D1)
            END IF
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            ND1=NIS1
            IF(ND1.GT.0) THEN
              CALL mma_allocate(D1,ND1,Label='D1')
              D1(:)=Zero
              CALL SPEC1D(IMLTOP,FACT,D2,D1)
            END IF
          END IF
          If (.NOT.ALLOCATED(D1)) CALL mma_allocate(D1,1,Label='D1')

          IF(ND1.GT.0) THEN
            XTST=DDOT_(ND1,D1,1,D1,1)
            IF(XTST.GT.1.0D12) THEN
              WRITE(6,'(1x,a,6i10)')' SIGMA G2 ICASE1,ISYM1:',
     &                                         ICASE1,ISYM1
              WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
     &                                         ICASE2,ISYM2
              Call Crash()
            END IF
          END IF

          DO ICASE2=ICASE1+1,NCASES
            IF(IFCOUP(ICASE2,ICASE1).EQ.0) Cycle
            DO ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) Cycle
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NSGMX=NIS2*NAS2
              IF(NSGMX.EQ.0) Cycle

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
                CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
                LSGMX=lg_SGMX
              ELSE
                LSGMX=Allocate_GA_Array(NSGMX,'SGMX')
              END IF

* SVC: this array is just zero....
*             XTST=DDOT_(NSGMX,GA_Array(LSGMX)%A,1,
*    &                         GA_Array(LSGMX)%A,1)
*             IF(XTST.GT.1.0D12) THEN
*               WRITE(6,'(1x,a,6i10)')' SIGMA H. ICASE1,ISYM1:',
*    &                                           ICASE1,ISYM1
*               WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
*    &                                           ICASE2,ISYM2
*               Call Crash()
*             END IF

#ifdef _DEBUGPRINT_
              WRITE(6,*)' ISYM1,ICASE1:',ISYM1,ICASE1
              WRITE(6,*)' ISYM2,ICASE2:',ISYM2,ICASE2
              WRITE(6,*)' SIGMA calling SGM with IMLTOP=',IMLTOP
#endif
C Compute contribution SGMX <- D2, and SGMX <- D1  if any
              CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &                 D1,D2,LSGMX,LISTS)

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                XTST=RHS_DDOT(NAS2,NIS2,lg_SGMX,lg_SGMX)
              ELSE
                XTST=DDOT_(NSGMX,GA_Arrays(LSGMX)%A,1,
     &                           GA_Arrays(LSGMX)%A,1)
              END IF

              IF(XTST.GT.1.0D12) THEN
                WRITE(6,'(1x,a,6i10)')' SIGMA I. ICASE1,ISYM1:',
     &                                           ICASE1,ISYM1
                WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
     &                                           ICASE2,ISYM2
                Call Crash()
              END IF

              IF (ICASE2.NE.12 .AND. ICASE2.NE.13) THEN
                MAX_MESG_SIZE = 2**27
                DO LSGMX_STA=1,NSGMX,MAX_MESG_SIZE
                  NSGMX_BLK=MIN(MAX_MESG_SIZE,NSGMX-LSGMX_STA+1)
                  CALL GADSUM(GA_Arrays(LSGMX)%A(LSGMX_STA),
     &                        NSGMX_BLK)
                END DO
                CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
                CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
                CALL RHS_ADD(NAS2,NIS2,lg_SGMX,GA_Arrays(LSGMX)%A)
                Call Deallocate_GA_Array(LSGMX)
              END IF

C-SVC: no need for the replicate arrays any more, fall back to one array
              CALL RHS_SAVE (NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
              CALL RHS_FREE (lg_SGMX)
            End Do
          End Do
          CALL mma_deallocate(D2)
          CALL mma_deallocate(D1)
        End Do
      End Do

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSGM=CPUSGM+(CPU1-CPU0)
      TIOSGM=TIOSGM+(TIO1-TIO0)

#ifdef _DEBUGPRINT_
      WRITE(6,*)' End of SIGMA. Flop counts:'
      WRITE(6,'(a,i12)')' In MLTSCA:',NFSCA
      WRITE(6,'(a,i12)')' In MLTDXP:',NFDXP
      WRITE(6,'(a,i12)')' In MLTMV :',NFMV
      WRITE(6,'(a,i12)')' In MLTR1 :',NFR1
      WRITE(6,*)
#endif

      Call mma_deallocate(FIT_Full)
      Call mma_deallocate(FTI_Full)
      Call mma_deallocate(FIA_Full)
      Call mma_deallocate(FAI_Full)
      Call mma_deallocate(FTA_Full)
      Call mma_deallocate(FAT_Full)
      Do iSym = 1, nSym
         nullify(FIT(iSym)%A,FTI(iSym)%A,FIA(iSym)%A,FAI(iSym)%A,
     &           FTA(iSym)%A,FAT(iSym)%A)
      End Do

C Transform contrav C  to eigenbasis of H0(diag):
      CALL PTRTOSR(1,IVEC,IVEC)
C Transform covar. sigma to eigenbasis of H0(diag):
      CALL PTRTOSR(0,JVEC,JVEC)

      CONTAINS
      Subroutine Crash()
         WRITE(6,*)' Colossal value detected in SIGMA.'
         WRITE(6,*)' This implies that the thresholds used for linear'
         WRITE(6,*)' dependence removal must be increased.'
         WRITE(6,*)' Present values, THRSHN, THRSHS:',THRSHN,THRSHS
         WRITE(6,*)' Use keyword THRESHOLD in input to increase these'
         WRITE(6,*)' values and then run again.'
         CALL ABEND()
      END Subroutine Crash
      END SUBROUTINE SIGMA_CASPT2
