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
* Copyright (C) 1988,1991,1992,1998, Per Ake Malmqvist                 *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE POLY3(IFF)
      use fciqmc_interface, only: DoFCIQMC
      use caspt2_global, only:iPrGlb
      use caspt2_global, only:LUCIEX, IDTCEX, LUSOLV
      use PrintLevel, only: VERBOSE
      use gugx, only: SGS, L2ACT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: DoCumulant, iSCF, jState, nActel,
     &                         nConf, nState, STSym, EPSA, mState
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_ || defined _DMRG_
      use caspt2_module, only: DMRG
#endif
      use pt2_guga, only: CIThr, nG1, nG2, nG3, nG3Tot, Eta
      IMPLICIT NONE
C  IBM TEST VERSION 0, 1988-06-23.
C  NEW VERSION 1991-02-23, FOR USE WITH RASSCF IN MOLCAS PACKAGE.
C  NEW VERSION 1992-12-05, FOR MOLCAS-3 VERSION.
C  NEW VERSION 1998-10-02
C  AUTHOR: PER-AAKE MALMQUIST

C THIS PROGRAM CALCULATES 1-EL, 2-EL, AND 3-EL
C DENSITY MATRICES FOR A CASSCF WAVE FUNCTION.
C IF THE INTEGER KEY IFF.EQ.1, THEN
C IT ALSO PRODUCES THE CONTRACTIONS OF 1-EL -- 4-EL
C DENSITY MATRICES WITH THE FOCK OPERATOR USED IN
C THE CASSCF-MP2 PROGRAM. THE RESULTS ARE WRITTEN
C TO FILE IN SEVERAL FORMS, TO SUPPORT BOTH KERSTINS
C PRESENT PROGRAM AND ALSO SUCH NEW PROCEDURES WHICH
C MIGHT TAKE ADVANTAGE OF ALL INDEX PERMUTATION SYMMETRIES.
C THE RDSTAT AND THE GUGA ROUTINES USED IN THIS
C PROGRAM ASSUMES THE JOBIPH IS PRODUCED BY THE RASSCF PROGRAM.


      INTEGER IFF

      INTEGER ILEV
      INTEGER NG3MAX
      INTEGER ILUID

      INTEGER IDCI
      INTEGER J

      INTEGER IPARDIV
      INTEGER*1, ALLOCATABLE :: idxG3(:,:)
      REAL*8, ALLOCATABLE, TARGET:: G1(:), G2(:), G3(:)
      REAL*8, ALLOCATABLE, TARGET:: F1_H(:), F2_H(:), F3_H(:)
      REAL*8, POINTER:: F1(:), F2(:), F3(:)
      REAL*8, ALLOCATABLE:: CI(:)

      Integer :: nLev
      nLev = SGS%nLev


      IF (IFF.EQ.1) THEN
C ORBITAL ENERGIES IN CI-COUPLING ORDER:
        DO ILEV=1,NLEV
          ETA(ILEV)=EPSA(L2ACT(ILEV))
        END DO
      END IF

      CALL mma_allocate(G1,NG1,LABEL='G1')
      CALL mma_allocate(G2,NG2,LABEL='G2')

C-SVC20100831: recompute approximate max NG3 size needed
      NG3MAX=iPARDIV(NG3TOT,NG2)

C-SVC20100831: allocate local G3 matrices
      CALL mma_allocate(G3,NG3MAX,LABEL='G3')

      CALL mma_allocate(idxG3,6,NG3MAX,label='idxG3')
      idxG3(:,:)=0

      G1(1)=0.0D0
      G2(1)=0.0D0
      G3(1)=0.0D0

C ALLOCATE SPACE FOR CORRESPONDING COMBINATIONS WITH H0:
      IF (IFF.EQ.1) THEN
        CALL mma_allocate(F1_H,NG1,LABEL='F1_H')
        CALL mma_allocate(F2_H,NG2,LABEL='F2_H')
        CALL mma_allocate(F3_H,NG3MAX,LABEL='F3_H')
        F1=>F1_H
        F2=>F2_H
        F3=>F3_H
      ELSE
        F1=>G1
        F2=>G2
        F3=>G3
      END IF

* NG3 will change inside subroutine MKFG3 to the actual
* number of nonzero elements, that is why here we allocate
* with NG3MAX, but we only store (PT2_PUT) the first NG3
* elements of the G3 and F3
      NG3=NG3MAX

      if (.not. DoFCIQMC) then
        CALL mma_allocate(CI,NCONF,Label='CI')

        IF (.NOT. DoCumulant .AND. ISCF.EQ.0) THEN
          IDCI=IDTCEX
          DO J=1,JSTATE-1
            CALL DDAFILE(LUCIEX,0,CI,NCONF,IDCI)
          END DO
          CALL DDAFILE(LUCIEX,2,CI,NCONF,IDCI)
          IF (IPRGLB.GE.VERBOSE) THEN
            WRITE(6,*)
            IF (NSTATE.GT.1) THEN
              WRITE(6,'(A,I4)')
     &       ' With new orbitals, the CI array of state ',MSTATE(JSTATE)
            ELSE
              WRITE(6,*)' With new orbitals, the CI array is:'
            END IF
            CALL PRWF_CP2(STSYM,NCONF,CI,CITHR)
          END IF
        ELSE
          CI(1)=1.0D0
        END IF
      end if

      IF (ISCF.NE.0.AND.NACTEL.NE.0) THEN
        CALL SPECIAL( G1,G2,G3,F1,F2,F3,idxG3)
      ELSE IF (ISCF.EQ.0) THEN
C-SVC20100903: during mkfg3, NG3 is set to the actual value
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_ || defined _DMRG_
        IF (.NOT. DoCumulant .AND. .NOT. DMRG) THEN
#endif
          If (.NOT.ALLOCATED(CI)) CALL mma_allocate(CI,1,LABEL='CI')
          CALL MKFG3(IFF,CI,G1,F1,G2,F2,G3,F3,idxG3,nLev)
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_ || defined _DMRG_
        ELSE
          CALL MKFG3DM(IFF,G1,F1,G2,F2,G3,F3,idxG3,nLev)
        END IF
#endif
      END IF

      if (ALLOCATED(CI)) CALL mma_deallocate(CI)

      IF(NLEV.GT.0) THEN
        CALL PT2_PUT(NG1,' GAMMA1',G1)
        CALL PT2_PUT(NG2,' GAMMA2',G2)
        CALL PT2_PUT(NG3,' GAMMA3',G3)
        iLUID=0
        CALL I1DAFILE(LUSOLV,1,idxG3,6*NG3,iLUID)
        IF(IFF.EQ.1) THEN
          CALL PT2_PUT(NG1,' DELTA1',F1)
          CALL PT2_PUT(NG2,' DELTA2',F2)
          CALL PT2_PUT(NG3,' DELTA3',F3)
        END IF
      END IF

      IF(NLEV.GT.0) THEN
        CALL mma_deallocate(G1)
        CALL mma_deallocate(G2)
        CALL mma_deallocate(G3)
        CALL mma_deallocate(idxG3)
        IF(IFF.EQ.1) THEN
          CALL mma_deallocate(F1_H)
          CALL mma_deallocate(F2_H)
          CALL mma_deallocate(F3_H)
        END IF
        F1=>Null()
        F2=>Null()
        F3=>Null()
      END IF

      END SUBROUTINE POLY3
