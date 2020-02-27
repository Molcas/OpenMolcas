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
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "sigma.fh"
#include "SysDef.fh"
      COMMON /CPLCAS/ IFCOUP(MXCASE,MXCASE)
      COMMON /FOCKOF/ LFIT,LFIA,LFTA,LFTI,LFAI,LFAT,
     &                IOFFIT(8),IOFFIA(8),IOFFTA(8)

C Compute |JVEC> := BETA* |JVEC> + ALPHA* (H0-E0)* |IVEC>
C where the vectors are represented in transformed basis and
C are  stored at positions IVEC and JVEC on the LUSOLV unit.

      CALL QENTER('SIGMA')

#ifdef _DEBUG_
      WRITE(6,*)' Entering SIGMA.'
      WRITE(6,*)
     &' Compute |JVEC> := Beta*|JVEC> + Alpha*(H0-E0)|IVEC>'
      WRITE(6,'(1x,a,2f15.6)')'Alpha,Beta:',Alpha,Beta
      WRITE(6,'(1x,a,2i5)')'IVEC,JVEC:',IVEC,JVEC
#endif

C Enter coupling cases for non-diagonal blocks:
      DO J=1,NCASES
      DO I=1,NCASES
      IFCOUP(I,J)=0
      END DO
      END DO
      IFCOUP( 2, 1)= 1
      IFCOUP( 3, 1)= 2
      IFCOUP( 5, 1)= 3
      IFCOUP( 6, 1)= 4
      IFCOUP( 7, 1)= 5
      IFCOUP( 6, 2)= 6
      IFCOUP( 7, 3)= 7
      IFCOUP( 5, 4)= 8
      IFCOUP( 8, 4)= 9
      IFCOUP( 9, 4)=10
      IFCOUP(10, 4)=11
      IFCOUP(11, 4)=12
      IFCOUP( 6, 5)=13
      IFCOUP( 7, 5)=14
      IFCOUP(10, 5)=15
      IFCOUP(11, 5)=16
      IFCOUP(12, 5)=23
      IFCOUP(13, 5)=24
      IFCOUP(12, 6)=17
      IFCOUP(13, 7)=18
      IFCOUP(10, 8)=19
      IFCOUP(11, 9)=20
      IFCOUP(12,10)=21
      IFCOUP(13,11)=22

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
      IF(ALPHA.EQ.0.0D0) GOTO 99
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
      NFIT=NFIT+1
      NFIA=NFIA+1
      NFTA=NFTA+1

      CALL GETMEM('FIT','ALLO','REAL',LFIT,NFIT)
      CALL GETMEM('FIA','ALLO','REAL',LFIA,NFIA)
      CALL GETMEM('FTA','ALLO','REAL',LFTA,NFTA)
      CALL GETMEM('FTI','ALLO','REAL',LFTI,NFIT)
      CALL GETMEM('FAI','ALLO','REAL',LFAI,NFIA)
      CALL GETMEM('FAT','ALLO','REAL',LFAT,NFTA)

      IFIFA=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NS=NSSH(ISYM)
        NO=NORB(ISYM)
        CALL FBLOCK(WORK(LFIFA+IFIFA),NO,NI,NA,NS,
     &       WORK(LFIT+IOFFIT(ISYM)),WORK(LFTI+IOFFIT(ISYM)),
     &       WORK(LFIA+IOFFIA(ISYM)),WORK(LFAI+IOFFIA(ISYM)),
     &       WORK(LFTA+IOFFTA(ISYM)),WORK(LFAT+IOFFTA(ISYM)))
        IFIFA=IFIFA+(NO*(NO+1))/2
      END DO

      CALL TIMING(CPU0,CPU,TIO0,TIO)
C Loop over types and symmetry block of sigma vector:
      DO 300 ICASE1=1,11
*     DO 300 ICASE1=1,NCASES
        DO 300 ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) GOTO 300
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          NSGM2=NIS1*NAS1
          IF(NSGM2.EQ.0) GOTO 300

          CALL GETMEM('SGM2','ALLO','REAL',LSGM2,NSGM2)
          CALL DCOPY_(NSGM2,[0.0D0],0,WORK(LSGM2),1)

          NSGM1=0
          LSGM1=1
          IF(ICASE1.EQ.1) THEN
            NSGM1=NASH(ISYM1)*NISH(ISYM1)
          ELSE IF(ICASE1.EQ.4) THEN
            NSGM1=NASH(ISYM1)*NSSH(ISYM1)
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            NSGM1=NIS1
          END IF
          IF(NSGM1.GT.0) THEN
            CALL GETMEM('SGM1','ALLO','REAL',LSGM1,NSGM1)
            CALL DCOPY_(NSGM1,[0.0D0],0,WORK(LSGM1),1)
          END IF

          IMLTOP=0
          DO 200 ICASE2=ICASE1+1,NCASES
            IFC=IFCOUP(ICASE2,ICASE1)
            IF(IFC.EQ.0) GOTO 200
            DO 100 ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) GOTO 100
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NCX=NIS2*NAS2
              IF(NCX.EQ.0) GOTO 100

              CALL RHS_ALLO(NAS2,NIS2,lg_CX)
              CALL RHS_READ(NAS2,NIS2,lg_CX,ICASE2,ISYM2,IVEC)
C SVC: for case H (12,13) we can now pass the distributed array ID to
C the SGM subroutines
              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                LCX=lg_CX
                XTST=RHS_DDOT(NAS2,NIS2,lg_CX,lg_CX)
              ELSE
                CALL GETMEM('CX','ALLO','REAL',LCX,NCX)
                CALL RHS_GET(NAS2,NIS2,lg_CX,WORK(LCX))
                CALL RHS_FREE(NAS2,NIS2,lg_CX)
                XTST=DDOT_(NCX,WORK(LCX),1,WORK(LCX),1)
              END IF

              IF(XTST.GT.1.0D12) THEN
                WRITE(6,'(1x,a,6i10)')' SIGMA A. ICASE2,ISYM2:',
     &                                           ICASE2,ISYM2
                GOTO 999
              END IF

#ifdef _DEBUG_
              WRITE(6,*)' ISYM1,ICASE1:',ISYM1,ICASE1
              WRITE(6,*)' ISYM2,ICASE2:',ISYM2,ICASE2
              WRITE(6,*)' SIGMA calling SGM with IMLTOP=',IMLTOP
#endif
C Compute contribution SGM2 <- CX, and SGM1 <- CX  if any
              CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &                 WORK(LSGM1),LSGM2,LCX,iWORK(LLISTS))

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                CALL RHS_FREE(NAS2,NIS2,lg_CX)
              ELSE
                CALL GETMEM('CX','FREE','REAL',LCX,NCX)
              END IF

C Check for colossal values of SGM2 and SGM1
              XTST=DDOT_(NSGM2,WORK(LSGM2),1,WORK(LSGM2),1)
              IF(XTST.GT.1.0D12) THEN
                WRITE(6,'(1x,a,6i10)')' SIGMA B. ICASE1,ISYM1:',
     &                                           ICASE1,ISYM1
                WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
     &                                           ICASE2,ISYM2
                GOTO 999
              END IF

              IF(NSGM1.GT.0) THEN
                XTST=DDOT_(NSGM1,WORK(LSGM1),1,WORK(LSGM1),1)
                IF(XTST.GT.1.0D12) THEN
                  WRITE(6,'(1x,a,6i10)')' SIGMA B2. ICASE1,ISYM1:',
     &                                              ICASE1,ISYM1
                  WRITE(6,'(1x,a,6i10)')'           ICASE2,ISYM2:',
     &                                              ICASE2,ISYM2
                  GOTO 999
                END IF
              END IF

 100        CONTINUE
 200      CONTINUE

C-SVC: sum the replicate arrays:
          MAX_MESG_SIZE = 2**27
          DO LSGM2_STA=1,NSGM2,MAX_MESG_SIZE
            NSGM2_BLK=MIN(MAX_MESG_SIZE,NSGM2-LSGM2_STA+1)
            CALL GADSUM(WORK(LSGM2+LSGM2_STA-1),NSGM2_BLK)
          END DO

          IF (NSGM1.GT.0) THEN
            CALL GADSUM(WORK(LSGM1),NSGM1)
          END IF

C       XTST2=DDOT_(NSGM2,WORK(LSGM2),1,WORK(LSGM2),1)
C       XTST1=0.0D0
C       IF(NSGM1.GT.0)XTST1=DDOT_(NSGM1,WORK(LSGM1),1,WORK(LSGM1),1)
C       WRITE(6,'(1x,a,a,i2,2f16.6)')
C    & 'Contr. SGM2, SGM1, ',cases(icase1),isym1,xtst2,xtst1

C If there are 1-electron contributions, add them into the 2-el
C part (This requires a non-empty active space.)
          IF(NSGM1.GT.0) THEN
            FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
            IF (ICASE1.EQ.1) THEN
              CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LSGM2),
     &                  WORK(LSGM1))
            ELSE IF(ICASE1.EQ.4) THEN
              CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LSGM2),
     &                  WORK(LSGM1))
            ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
              CALL SPEC1D(IMLTOP,FACT,WORK(LSGM2),WORK(LSGM1))
            END IF

            XTST=DDOT_(NSGM2,WORK(LSGM2),1,WORK(LSGM2),1)
            IF(XTST.GT.1.0D12) THEN
              WRITE(6,'(1x,a,6i10)')' SIGMA C. ICASE1,ISYM1:',
     &                                         ICASE1,ISYM1
              GOTO 999
            END IF

            CALL GETMEM('SGM1','FREE','REAL',LSGM1,NSGM1)
          END IF

C-SVC: no need for the replicate arrays any more, fall back to one array
          CALL RHS_ALLO (NAS1,NIS1,lg_SGM2)
          CALL RHS_PUT (NAS1,NIS1,lg_SGM2,WORK(LSGM2))
          CALL GETMEM('SGM2','FREE','REAL',LSGM2,NSGM2)

C Add to sigma array. Multiply by S to  lower index.
          NSGMX=NSGM2
          CALL RHS_ALLO(NAS1,NIS1,lg_SGMX)
          CALL RHS_READ(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)

          XTST=RHS_DDOT(NAS1,NIS1,lg_SGMX,lg_SGMX)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA D. ICASE1,ISYM1:',ICASE1,ISYM1
            WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',ICASE2,ISYM2
            GOTO 999
          END IF

*         IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
            CALL RHS_STRANS(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX,
     &                      ICASE1,ISYM1)
*         ELSE
*           CALL RHS_DAXPY(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX)
*         END IF
          CALL RHS_FREE (NAS1,NIS1,lg_SGM2)

          XTST=RHS_DDOT(NAS1,NIS1,lg_SGMX,lg_SGMX)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA E. ICASE1,ISYM1:',ICASE1,ISYM1
            WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',ICASE2,ISYM2
            GOTO 999
          END IF

C Write SGMX to disk.
          CALL RHS_SAVE (NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
          CALL RHS_FREE (NAS1,NIS1,lg_SGMX)
 300  CONTINUE

      IMLTOP=1
C Loop over types and symmetry block of CX vector:
      DO 600 ICASE1=1,11
*     DO 600 ICASE1=1,NCASES
        DO 600 ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) GOTO 600
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          ND2=NIS1*NAS1
          IF(ND2.EQ.0) GOTO 600

          CALL RHS_ALLO (NAS1,NIS1,lg_D2)
          CALL RHS_SCAL (NAS1,NIS1,lg_D2,0.0D0)
C Contract S*CX to form D2. Also form D1 from D2, if needed.

          NCX=ND2
          CALL RHS_ALLO (NAS1,NIS1,lg_CX)
          CALL RHS_READ (NAS1,NIS1,lg_CX,ICASE1,ISYM1,IVEC)

          XTST=RHS_DDOT(NAS1,NIS1,lg_CX,lg_CX)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA F. ICASE1,ISYM1:',ICASE1,ISYM1
            GOTO 999
          END IF

          IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
           CALL RHS_STRANS (NAS1,NIS1,ALPHA,lg_CX,lg_D2,ICASE1,ISYM1)
          ELSE
           CALL RHS_DAXPY(NAS1,NIS1,ALPHA,lg_CX,lg_D2)
          END IF
          CALL RHS_FREE (NAS1,NIS1,lg_CX)

CPAM Sanity check:
          XTST=RHS_DDOT(NAS1,NIS1,lg_D2,lg_D2)
          IF(XTST.GT.1.0D12) THEN
            WRITE(6,'(1x,a,6i10)')' SIGMA G1 ICASE1,ISYM1:',ICASE1,ISYM1
            WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',ICASE2,ISYM2
            GOTO 999
          END IF

          CALL GETMEM('D2','ALLO','REAL',LD2,ND2)
          CALL RHS_GET (NAS1,NIS1,lg_D2,WORK(LD2))
          CALL RHS_FREE (NAS1,NIS1,lg_D2)

          ND1=0
          LD1=1
          IMLTOP=1
          FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
          IF(ICASE1.EQ.1) THEN
            ND1=NASH(ISYM1)*NISH(ISYM1)
            IF(ND1.GT.0) THEN
              CALL GETMEM('D1','ALLO','REAL',LD1,ND1)
              CALL DCOPY_(ND1,[0.0D0],0,WORK(LD1),1)
              CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LD2),
     &                    WORK(LD1))
            END IF
          ELSE IF(ICASE1.EQ.4) THEN
            ND1=NASH(ISYM1)*NSSH(ISYM1)
            IF(ND1.GT.0) THEN
              CALL GETMEM('D1','ALLO','REAL',LD1,ND1)
              CALL DCOPY_(ND1,[0.0D0],0,WORK(LD1),1)
              CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LD2),
     &                    WORK(LD1))
            END IF
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            ND1=NIS1
            IF(ND1.GT.0) THEN
              CALL GETMEM('D1','ALLO','REAL',LD1,ND1)
              CALL DCOPY_(ND1,[0.0D0],0,WORK(LD1),1)
              CALL SPEC1D(IMLTOP,FACT,WORK(LD2),WORK(LD1))
            END IF
          END IF

          IF(ND1.GT.0) THEN
            XTST=DDOT_(ND1,WORK(LD1),1,WORK(LD1),1)
            IF(XTST.GT.1.0D12) THEN
              WRITE(6,'(1x,a,6i10)')' SIGMA G2 ICASE1,ISYM1:',
     &                                         ICASE1,ISYM1
              WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
     &                                         ICASE2,ISYM2
              GOTO 999
            END IF
          END IF

          DO 500 ICASE2=ICASE1+1,NCASES
            IF(IFCOUP(ICASE2,ICASE1).EQ.0) GOTO 500
            DO 400 ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) GOTO 400
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NSGMX=NIS2*NAS2
              IF(NSGMX.EQ.0) GOTO 400

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
                CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
                LSGMX=lg_SGMX
              ELSE
                CALL GETMEM('SGMX','ALLO','REAL',LSGMX,NSGMX)
                CALL DCOPY_(NSGMX,[0.0D0],0,WORK(LSGMX),1)
              END IF

* SVC: this array is just zero....
*             XTST=DDOT_(NSGMX,WORK(LSGMX),1,WORK(LSGMX),1)
*             IF(XTST.GT.1.0D12) THEN
*               WRITE(6,'(1x,a,6i10)')' SIGMA H. ICASE1,ISYM1:',
*    &                                           ICASE1,ISYM1
*               WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
*    &                                           ICASE2,ISYM2
*               GOTO 999
*             END IF

#ifdef _DEBUG_
              WRITE(6,*)' ISYM1,ICASE1:',ISYM1,ICASE1
              WRITE(6,*)' ISYM2,ICASE2:',ISYM2,ICASE2
              WRITE(6,*)' SIGMA calling SGM with IMLTOP=',IMLTOP
#endif
C Compute contribution SGMX <- D2, and SGMX <- D1  if any
              CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &                 WORK(LD1),LD2,LSGMX,iWORK(LLISTS))

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                XTST=RHS_DDOT(NAS2,NIS2,lg_SGMX,lg_SGMX)
              ELSE
                XTST=DDOT_(NSGMX,WORK(LSGMX),1,WORK(LSGMX),1)
              END IF

              IF(XTST.GT.1.0D12) THEN
                WRITE(6,'(1x,a,6i10)')' SIGMA I. ICASE1,ISYM1:',
     &                                           ICASE1,ISYM1
                WRITE(6,'(1x,a,6i10)')'          ICASE2,ISYM2:',
     &                                           ICASE2,ISYM2
                GOTO 999
              END IF

              IF (ICASE2.NE.12 .AND. ICASE2.NE.13) THEN
                MAX_MESG_SIZE = 2**27
                DO LSGMX_STA=1,NSGMX,MAX_MESG_SIZE
                  NSGMX_BLK=MIN(MAX_MESG_SIZE,NSGMX-LSGMX_STA+1)
                  CALL GADSUM(WORK(LSGMX+LSGMX_STA-1),NSGMX_BLK)
                END DO
                CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
                CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
                CALL RHS_ADD(NAS2,NIS2,lg_SGMX,WORK(LSGMX))
                CALL GETMEM('SGMX','FREE','REAL',LSGMX,NSGMX)
              END IF

C-SVC: no need for the replicate arrays any more, fall back to one array
              CALL RHS_SAVE (NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
              CALL RHS_FREE (NAS2,NIS2,lg_SGMX)
 400        CONTINUE
 500      CONTINUE
          CALL GETMEM('D2','FREE','REAL',LD2,ND2)
          IF(ND1.GT.0) CALL GETMEM('D1','FREE','REAL',LD1,ND1)
 600  CONTINUE

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSGM=CPUSGM+(CPU1-CPU0)
      TIOSGM=TIOSGM+(TIO1-TIO0)

#ifdef _DEBUG_
      WRITE(6,*)' End of SIGMA. Flop counts:'
      WRITE(6,'(a,i12)')' In MLTSCA:',NFSCA
      WRITE(6,'(a,i12)')' In MLTDXP:',NFDXP
      WRITE(6,'(a,i12)')' In MLTMV :',NFMV
      WRITE(6,'(a,i12)')' In MLTR1 :',NFR1
      WRITE(6,*)
#endif

      CALL GETMEM('FIT','FREE','REAL',LFIT,NFIT)
      CALL GETMEM('FIA','FREE','REAL',LFIA,NFIA)
      CALL GETMEM('FTA','FREE','REAL',LFTA,NFTA)
      CALL GETMEM('FTI','FREE','REAL',LFTI,NFIT)
      CALL GETMEM('FAI','FREE','REAL',LFAI,NFIA)
      CALL GETMEM('FAT','FREE','REAL',LFAT,NFTA)

C Transform contrav C  to eigenbasis of H0(diag):
      CALL PTRTOSR(1,IVEC,IVEC)
C Transform covar. sigma to eigenbasis of H0(diag):
      CALL PTRTOSR(0,JVEC,JVEC)

  99  CONTINUE
      CALL QEXIT('SIGMA')
      RETURN

 999  CONTINUE
C Error exit.
      WRITE(6,*)' Colossal value detected in SIGMA.'
      WRITE(6,*)' This implies that the thresholds used for linear'
      WRITE(6,*)' dependence removal must be increased.'
      WRITE(6,*)' Present values, THRSHN, THRSHS:',THRSHN,THRSHS
      WRITE(6,*)' Use keyword THRESHOLD in input to increase these'
      WRITE(6,*)' values and then run again.'
      CALL ABEND()
      END
