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
      SUBROUTINE SIGDER(IVEC,JVEC,SCAL)
      use definitions, only: wp, iwp
      use constants, only: Zero, One, Two
      use Fockof, only: FAI_Full, FAT_Full, FIA_Full, FIT_Full,
     &                  FTA_Full, FTI_Full, IOFFIT, IOFFIA, IOFFTA,
     &                  FIT, FTI, FIA, FAI, FAT, FTA
      use caspt2_global, only: LUSTD,idSDMat
      use caspt2_global, only: FIFA,LISTS
      use stdalloc, only: mma_allocate, mma_deallocate
      use EQSOLV, only: IFCOUP
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array,
     &                   GA_Arrays
      use caspt2_module, only: CPUSGM, TIOSGM, FockType, G1SecIn,
     &                         nActEl, nCases, nSym, nASup, nIsh, nAsh,
     &                         nSsh, nOrb, nInDep, nISup
      IMPLICIT None
#if defined(_MOLCAS_MPP_) && defined(_GA_)
#include "global.fh"
#endif
      integer(kind=iwp), intent(in):: iVec, jVec
      real(kind=wp), intent(in) :: Scal

      real(kind=wp) CPU, CPU0, CPU1
      real(kind=wp) TIO, TIO0, TIO1
      real(kind=wp) FACT
      integer(kind=iwp) iCase, iCase1, idSDer, IfC, IFIFA, iLoop,
     &                  IMLTOP, ISYM, ISYM1, ISYM2, lCX, lg_CX, lg_D2,
     &                  lg_Sgm2, Lg_SgmX, lg_V1, lSgmX,
     &                  lSgmX_Sta, Max_MESG_SIZE, NA, NAS, NAS1, NAS2,
     &                  nCX, nD1, NFIA, NFIT, NFTA, NI, NIN1, NIN2,
     &                  NIS1, NIS2, nLoop, NO, NS, nSgm1, nSgm2,
     &                  lSgm2_Sta, nD2, iCase2, nSgm2_Blk, nSgmX,
     &                  nSgmX_Blk
      real(kind=wp), allocatable :: WRK(:),SDER1(:),SDER2(:),SGM1(:),
     &                              D1(:), SGM2(:), D2(:)
      integer(kind=iwp) MaxLen
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      logical(kind=iwp) :: bStat
#endif
C
C     Work in the MO basis
C     We need both explicit and implicit overlap derivatives. The latter
C     comes from the derivative of the transformation matrix.
C
C     p,q: inactive or secondary
C     y,z: active (t,u)
C     a,b: internally contracted
C     T1_{px} * S1_{xy} f_{yz} * T2_{pz}
C     = T1pa*C1xa * S1xy * fyz * T2pb*C2zb
C     Derivative of S1:
C     = (T1Ct1)px * (T2Ct2*f)py * dS1xy/da
C     Derivative of C1 (or, Lagrangian multiplier in MO basis):
C     = T1pa*dC1xa/da * S1xy * fyz * T2pb*C2zb
C       ...
C     = -1/2 (T1Ct1)pu * dS1tu/da * (T2Ct2*f*S1C1*Ct1)pt
C     Derivative of C2 (or, Lagrangian multiplier in MO basis):
C     = T1pa*C1xa * S1xy * fyz * T2pb*dC2zb/da
C       ...
C     = -1/2 (T1Ct1St1*f*C2*Ct2)pt * (T2Ct2)pu * dS2tu/da
C
C     About IMLTOP for the SGM subroutine
C     With IMLTOP=0: the vector for the second argument has to be
C     contravariant form (T*C),
C     With IMLTOP=1: the vector for the first  argument has to be
C     covariant form (T*SC),
C
C
C     Allocate some matrices for storing overlap and transformation
C     derivatives. Here constructs these derivatives in the MO basis,
C     but not in the internally contracted basis.
C
      MaxLen = 0
      Do iCase = 1, 11
        Do iSym = 1, nSym
          nAS = nASUP(iSym,iCase)
          MaxLen = Max(MaxLen,nAS*nAS)
        End Do
      End Do
C
      call mma_allocate(WRK,MaxLen,Label='WRK')
      Call DCopy_(MaxLen,[Zero],0,WRK,1)
C
      Do iCase = 1, 11
        Do iSym = 1, nSym
          nAS = nASUP(iSym,iCase)
          idSDer = idSDMat(iSym,iCase)
          CALL DDAFILE(LuSTD,1,WRK,nAS*nAS,idSDer)
        End Do
      End Do
      call mma_deallocate(WRK)
C
C If the G1 correction to the Fock matrix is used, then the
C inactive/virtual coupling elements (which are non-zero for the
C case of average CASSCF) cannot be used in the CASPT2 equations.
      IF(FOCKTYPE.EQ.'G1      ' .AND. (.NOT. G1SECIN)) THEN
        IFCOUP(12,5)=0
        IFCOUP(13,5)=0
      END IF


C Transform to standard representation:
      CALL PTRTOC(0,IVEC,IVEC) !! T*C (internally contracted -> MO)
      IF(IVEC.NE.JVEC) CALL PTRTOC(0,JVEC,JVEC)

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
      NFIT=NFIT+1
      NFIA=NFIA+1
      NFTA=NFTA+1

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

      END DO

      CALL TIMING(CPU0,CPU,TIO0,TIO)
C
C     Is it possible to reduce to one loop? We have to compute bra and
C     ket overlap and bra and ket wavefunctions are not identical, so
C     it seems impossible to reduce?
C
      NLOOP=2
      DO ILOOP=1,NLOOP
        !! ILOOP1 : <T+lambda|H|T       >
        !! ILOOP2 : <T       |H|T+lambda>

C Loop over types and symmetry block of sigma vector:
      DO ICASE1=1,11
*     DO ICASE1=1,NCASES
        DO ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) Cycle
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          NIN1=NINDEP(ISYM1,ICASE1)
          NSGM2=NIS1*NAS1
          IF(NSGM2.EQ.0) Cycle

          CALL mma_allocate(SGM2,NSGM2,Label='SGM2')
          SGM2(:)=Zero

          NSGM1=0
C         LSGM1=1
          IF(ICASE1.EQ.1) THEN
            NSGM1=NASH(ISYM1)*NISH(ISYM1)
          ELSE IF(ICASE1.EQ.4) THEN
            NSGM1=NASH(ISYM1)*NSSH(ISYM1)
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            NSGM1=NIS1
          END IF
          CALL mma_allocate(SGM1,MAX(1,NSGM1),Label='SGM1')
          SGM1(:) = Zero

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
              IF (IVEC.NE.JVEC .AND. ILOOP.EQ.2) THEN
                !! T = T + \lambda
                If (SCAL.ne.One)
     &             CALL RHS_SCAL(NAS2,NIS2,lg_CX,SCAL)
                CALL RHS_ALLO(NAS2,NIS2,lg_V1)
                CALL RHS_READ(NAS2,NIS2,lg_V1,ICASE2,ISYM2,JVEC)
                CALL RHS_DAXPY(NAS2,NIS2,One,lg_V1,lg_CX)
                CALL RHS_FREE(lg_V1)
              END IF
C SVC: for case H (12,13) we can now pass the distributed array ID to
C the SGM subroutines
              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                LCX=lg_CX
C               XTST=RHS_DDOT(NAS2,NIS2,lg_CX,lg_CX)
              ELSE
                LCX=Allocate_GA_Array(NCX,'CX')
                CALL RHS_GET(NAS2,NIS2,lg_CX,GA_Arrays(LCX)%A)
                CALL RHS_FREE(lg_CX)
C               XTST=DDOT_(NCX,GA_Arrays(LCX)%A,1,
C    &                         GA_Arrays(LCX)%A,1)
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
            End Do
          End Do

C-SVC: sum the replicate arrays:
          MAX_MESG_SIZE = 2**27
          DO LSGM2_STA=1,NSGM2,MAX_MESG_SIZE
            NSGM2_BLK=MIN(MAX_MESG_SIZE,NSGM2-LSGM2_STA+1)
            CALL GADSUM(SGM2(LSGM2_STA),NSGM2_BLK)
          END DO

          IF (NSGM1.GT.0) THEN
            CALL GADSUM(SGM1,NSGM1)
          END IF

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
          END IF
          call mma_deallocate(SGM1)

C-SVC: no need for the replicate arrays any more, fall back to one array
          CALL RHS_ALLO (NAS1,NIS1,lg_SGM2)
          CALL RHS_PUT (NAS1,NIS1,lg_SGM2,SGM2)
          CALL mma_deallocate(SGM2)

C Add to sigma array. Multiply by S to  lower index.
C         CALL RHS_ALLO(NAS1,NIS1,lg_SGMX)
C         CALL RHS_READ(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
          IF (ICASE1.LE.11) THEN
            CALL RHS_ALLO(NAS1,NIS1,lg_CX)
            CALL RHS_READ(NAS1,NIS1,lg_CX,ICASE1,ISYM1,IVEC)
            If (IVEC.NE.JVEC .AND. ILOOP.EQ.1) Then
              !! T = T + \lambda
              If (SCAL.ne.One) CALL RHS_SCAL(NAS1,NIS1,lg_CX,SCAL)
              CALL RHS_ALLO(NAS1,NIS1,lg_V1)
              CALL RHS_READ(NAS1,NIS1,lg_V1,ICASE1,ISYM1,JVEC)
              CALL RHS_DAXPY(NAS1,NIS1,One,lg_V1,lg_CX)
              CALL RHS_FREE(lg_V1)
            End If

            call mma_allocate(SDER1,NAS1*NAS1,Label='SDER1')
            idSDer = idSDMat(iSym1,iCase1)
            CALL DDAFILE(LuSTD,2,SDER1,nAS1*nAS1,idSDer)

            Call C1S1DER(SDER1)

            idSDer = idSDMat(iSym1,iCase1)
            CALL DDAFILE(LuSTD,1,SDER1,nAS1*nAS1,idSDer)
            call mma_deallocate(SDER1)

            CALL RHS_FREE(lg_CX)
          END IF

*         IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
C           CALL RHS_STRANS(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX,
C    &                      ICASE1,ISYM1)
*         ELSE
*           CALL RHS_DAXPY(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX)
*         END IF
          CALL RHS_FREE (lg_SGM2)

C Write SGMX to disk.
C         CALL RHS_SAVE (NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
C         CALL RHS_FREE (lg_SGMX)
        End Do
      End Do

      IMLTOP=1
C Loop over types and symmetry block of CX vector:
      DO ICASE1=1,11
*     DO ICASE1=1,NCASES
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

          IF (IVEC.NE.JVEC .AND. ILOOP.EQ.1) THEN
            !! T = T + \lambda
            If (SCAL.ne.One) CALL RHS_SCAL(NAS1,NIS1,lg_CX,SCAL)
            CALL RHS_ALLO(NAS1,NIS1,lg_V1)
            CALL RHS_READ(NAS1,NIS1,lg_V1,ICASE1,ISYM1,JVEC)
            CALL RHS_DAXPY(NAS1,NIS1,One,lg_V1,lg_CX)
            CALL RHS_FREE(lg_V1)
          END IF

          IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
           CALL RHS_STRANS(NAS1,NIS1,One,lg_CX,lg_D2,
     &                     ICASE1,ISYM1)
          ELSE
           CALL RHS_DAXPY(NAS1,NIS1,One,lg_CX,lg_D2)
          END IF
          CALL RHS_FREE (lg_CX)

          CALL mma_allocate(D2,ND2,Label='D2')
          CALL RHS_GET (NAS1,NIS1,lg_D2,D2)
          CALL RHS_FREE (lg_D2)

          ND1=0
C         LD1=1
          IMLTOP=1
          FACT=One/(DBLE(MAX(1,NACTEL)))
          IF(ICASE1.EQ.1) THEN
            ND1=NASH(ISYM1)*NISH(ISYM1)
            IF(ND1.GT.0) THEN
              call mma_allocate(D1,ND1,Label='D1')
              D1(:) = Zero
              CALL SPEC1A(IMLTOP,FACT,ISYM1,D2,D1)
            END IF
          ELSE IF(ICASE1.EQ.4) THEN
            ND1=NASH(ISYM1)*NSSH(ISYM1)
            IF(ND1.GT.0) THEN
              call mma_allocate(D1,ND1,Label='D1')
              D1(:) = Zero
              CALL SPEC1C(IMLTOP,FACT,ISYM1,D2,D1)
            END IF
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            ND1=NIS1
            IF(ND1.GT.0) THEN
              call mma_allocate(D1,ND1,Label='D1')
              D1(:) = Zero
              CALL SPEC1D(IMLTOP,FACT,D2,D1)
            END IF
          END IF
          If (.NOT.ALLOCATED(D1)) CALL mma_allocate(D1,1,Label='D1')

          !! No need to compute for ICASE2 = 12 and 13
          DO ICASE2=ICASE1+1,11 !! NCASES
            IF(IFCOUP(ICASE2,ICASE1).EQ.0) Cycle
            DO ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) Cycle
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NIN2=NINDEP(ISYM2,ICASE2)
              NSGMX=NIS2*NAS2
              IF(NSGMX.EQ.0) Cycle

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
                CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
                LSGMX=lg_SGMX
              ELSE
                LSGMX=Allocate_GA_Array(NSGMX,'SGMX')
              END IF

#ifdef _DEBUGPRINT_
              WRITE(6,*)' ISYM1,ICASE1:',ISYM1,ICASE1
              WRITE(6,*)' ISYM2,ICASE2:',ISYM2,ICASE2
              WRITE(6,*)' SIGMA calling SGM with IMLTOP=',IMLTOP
#endif
C Compute contribution SGMX <- D2, and SGMX <- D1  if any
              CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &                 D1,D2,LSGMX,LISTS)
C             If (iCase2.LE.11) Then
C             End If

              IF (ICASE2.NE.12 .AND. ICASE2.NE.13) THEN
                MAX_MESG_SIZE = 2**27
                DO LSGMX_STA=1,NSGMX,MAX_MESG_SIZE
                  NSGMX_BLK=MIN(MAX_MESG_SIZE,NSGMX-LSGMX_STA+1)
                  CALL GADSUM(GA_Arrays(LSGMX)%A(LSGMX_STA),
     &                        NSGMX_BLK)
                END DO
C               CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
C               CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
C               CALL RHS_ADD(NAS2,NIS2,lg_SGMX,GA_Array(LSGMX)%A)
                !! do C2DER
                call mma_allocate(SDER2,NAS2*NAS2,Label='SDER2')
                idSDer = idSDMat(iSym2,iCase2)
                CALL DDAFILE(LuSTD,2,SDER2,nAS2*nAS2,idSDer)

                Call C2DER(SDER2)

                idSDer = idSDMat(iSym2,iCase2)
                CALL DDAFILE(LuSTD,1,SDER2,nAS2*nAS2,idSDer)
                call mma_deallocate(SDER2)
                !!
C               Call Deallocate_GA_Array(LSGMX)
              END IF

C-SVC: no need for the replicate arrays any more, fall back to one array
C             CALL RHS_SAVE (NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
              IF (ICASE2.EQ.12 .OR.ICASE2.EQ.13) THEN
                CALL RHS_FREE (lg_SGMX)
              ELSE
                Call Deallocate_GA_Array(LSGMX)
              END IF
            End Do
          End Do
          CALL mma_deallocate(D2)
          call mma_deallocate(D1)
        End Do
      End Do

      End Do
C
C
C
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSGM=CPUSGM+(CPU1-CPU0)
      TIOSGM=TIOSGM+(TIO1-TIO0)
C
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

C Transform contrav C  to ei%Agenbasis of H0(diag):
      CALL PTRTOSR(1,IVEC,IVEC)
      IF(IVEC.NE.JVEC) CALL PTRTOSR(1,JVEC,JVEC)


      Contains
C
C-----------------------------------------------------------------------
C
      Subroutine C1S1DER(SDER)
C
      Implicit None
C
      real(kind=wp), intent(inout):: SDER(*)

      integer(kind=iwp) iType
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      integer(kind=iwp) lg_SDER
#endif
C
C     (T2Ct2*f)py * (T1Ct1)pz * dS1yz/da
C     -1/2 (T2Ct2*f*S1*C1*Ct1)pt * (T1Ct1)pu * dS1tu/da
C
      !! Finalize the derivative of S1
      !! 2S. (T2Ct2*f) * T1Ct1
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        CALL GA_CREATE_STRIPED ('H',NAS1,NAS1,'SDER',lg_SDER)
        CALL GA_PUT(lg_SDER,1,NAS1,1,NAS1,SDER,NAS1)
        call GA_DGEMM ('N','T',NAS1,NAS1,NIS1,
     *                 Two,lg_CX,lg_SGM2,One,lg_SDER)
      else
#endif
        Call DGEMM_('N','T',NAS1,NAS1,NIS1,
     *              Two,GA_Arrays(lg_CX)%A,NAS1,
     &                      GA_Arrays(lg_SGM2)%A,NAS1,
     *              One,SDER,NAS1)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
C     do i = 1, nas1*nis1
C       write (*,'(i4,2f20.10)') ,i,GA_Arrays(lg_cx)%A(i),
C    &                              GA_Arrays(lg_sgm2)%A(i)
C     end do
C
      !! Next, the derivative of C1
      !! 2C. (T2Ct2*f) * S1*C1 (MO -> IC)
      !!     lg_T * lg_V2 -> lg_V1
      CALL RHS_ALLO(NIN1,NIS1,lg_V1)
      ITYPE=1
      !! SGM2 is local quantity, so put this in GA?
      CALL RHS_SR2C (ITYPE,1,NAS1,NIS1,NIN1,lg_V1,lg_SGM2,
     &               ICASE1,ISYM1)
      !! 3C. (T2Ct2*f) * S1*C1 * Ct1 (IC -> MO)
      !!     lg_T * lg_V1 -> lg_V2
      ITYPE=0
      CALL RHS_SR2C (ITYPE,0,NAS1,NIS1,NIN1,lg_V1,lg_SGM2,
     &               ICASE1,ISYM1)
      CALL RHS_FREE(lg_V1)
C
      !! 4C. (T1Ct1*f) * (T2Ct2St2*f*C1*Ct1)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        call GA_DGEMM ('N','T',NAS1,NAS1,NIS1,
     *                -One,lg_CX,lg_SGM2,One,lg_SDER)
        CALL GA_GET(lg_SDER,1,NAS1,1,NAS1,SDER,NAS1)
        bStat = GA_destroy(lg_SDER)
      else
#endif
        Call DGEMM_('N','T',NAS1,NAS1,NIS1,
     *             -One,GA_Arrays(lg_CX)%A,NAS1,
     &                      GA_Arrays(lg_SGM2)%A,NAS1,
     *              One,SDER,NAS1)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
C
      End Subroutine C1S1DER
C
C-----------------------------------------------------------------------
C
      Subroutine C2DER(SDER)
C
      Implicit None
C
      real(kind=wp), intent(inout):: SDER(*)

      integer(kind=iwp) iType, lg_Sgm, lg_v2
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      integer(kind=iwp) lg_SDER
#endif
C
C     -1/2 (T2Ct2)pu * dS2tu/da * (T1Ct1St1*f*C2*Ct2)pt
C
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        CALL GA_CREATE_STRIPED ('V',NAS2,NIS2,'SDER',lg_SGMX)
        CALL GA_PUT(lg_SGMX,1,NAS2,1,NIS2,GA_Arrays(LSGMX)%A,NAS2)
      else
#endif
       lg_SGMX = LSGMX
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
C
      !! For icase = 12 or 13, there is no need to transform,
      !! so LTMP is always replicated, but T for A and C are
      !! distributed?, so...
      CALL RHS_ALLO(NIN2,NIS2,lg_V2)
      !! 2. (T1Ct1St1*f) * C2 (MO -> IC; LTMP -> LTMP2)
      ITYPE=0
      CALL RHS_SR2C (ITYPE,1,NAS2,NIS2,NIN2,lg_V2,lg_SGMX,
     &               ICASE2,ISYM2)
      !! 3. (T1Ct1St1*f) * C2 * Ct2 (IC -> MO; LTMP2 -> LTMP)
      CALL RHS_SR2C (ITYPE,0,NAS2,NIS2,NIN2,lg_V2,lg_SGMX,
     &               ICASE2,ISYM2)
      CALL RHS_FREE(lg_V2)
C
      !! 4. (T2Ct2*f) * (T1Ct1St1*f*C2*Ct2)
      CALL RHS_ALLO(NAS2,NIS2,lg_SGM)
      CALL RHS_READ(NAS2,NIS2,lg_SGM,ICASE2,ISYM2,IVEC)
          IF (IVEC.NE.JVEC .AND. ILOOP.EQ.2) THEN
            !! T = T + \lambda
            If (SCAL.ne.One) CALL RHS_SCAL(NAS2,NIS2,lg_SGM,SCAL)
            CALL RHS_ALLO(NAS2,NIS2,lg_V1)
            CALL RHS_READ(NAS2,NIS2,lg_V1,ICASE2,ISYM2,JVEC)
            CALL RHS_DAXPY(NAS2,NIS2,One,lg_V1,lg_SGM)
            CALL RHS_FREE(lg_V1)
          END IF
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        CALL GA_CREATE_STRIPED ('H',NAS2,NAS2,'SDER',lg_SDER)
        CALL GA_PUT(lg_SDER,1,NAS2,1,NAS2,SDER,NAS2)
        call GA_DGEMM ('N','T',NAS2,NAS2,NIS2,
     *                -One,lg_SGM,lg_SGMX,One,lg_SDER)
        CALL GA_GET(lg_SDER,1,NAS2,1,NAS2,SDER,NAS2)
        bStat = GA_destroy(lg_SGMX)
        bStat = GA_destroy(lg_SDER)
      else
#endif
        Call DGEMM_('N','T',NAS2,NAS2,NIS2,
     *             -One,GA_Arrays(lg_SGM)%A,NAS2,
     &                      GA_Arrays(LSGMX)%A,NAS2,
     *              One,SDER,NAS2)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
      CALL RHS_FREE(lg_SGM)
C
      End Subroutine C2DER
C
      End subroutine sigder
