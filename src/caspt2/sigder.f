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
C     SUBROUTINE TRDNS2O(IVEC,JVEC,DPT2)
      SUBROUTINE SIGDER(IVEC,JVEC,SCAL)
      use Fockof
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "sigma.fh"
#include "SysDef.fh"
#include "caspt2_grad.fh"
      COMMON /CPLCAS/ IFCOUP(MXCASE,MXCASE)
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
      Call GETMEM('WRK','ALLO','REAL',ipWRK,MaxLen)
      Call DCopy_(MaxLen,[0.0D+00],0,Work(ipWRK),1)
C
      idSD = 1
      Do iCase = 1, 11
        Do iSym = 1, nSym
          idSDMat(iSym,iCase) = idSD
          nAS = nASUP(iSym,iCase)
          CALL DDAFILE(LuSTD,0,Work(ipWRK),nAS*nAS,idSD)
          idSDer = idSDMat(iSym,iCase)
          ! idSDMat(iSym,iCase))
          CALL DDAFILE(LuSTD,1,Work(ipWRK),nAS*nAS,idSDer)
        End Do
      End Do
      Call GETMEM('WRK','FREE','REAL',ipWRK,MaxLen)
C
      Call GETMEM('SDER1','ALLO','REAL',ipSDER1,MaxLen)
      Call GETMEM('SDER2','ALLO','REAL',ipSDER2,MaxLen)
C
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

      IFIFA=0
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

        CALL FBLOCK(WORK(LFIFA+IFIFA),NO,NI,NA,NS,
     &              FIT(ISYM)%A(:),FTI(ISYM)%A(:),
     &              FIA(ISYM)%A(:),FAI(ISYM)%A(:),
     &              FTA(ISYM)%A(:),FAT(ISYM)%A(:))

        IFIFA=IFIFA+(NO*(NO+1))/2

      END DO

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      NLOOP=2
      DO 1000 ILOOP=1,NLOOP
        ! IF(ILOOP.EQ.1) THEN
        !   IBRA=IVEC
        !   IKET=JVEC
        ! ELSE
        !   IBRA=JVEC
        !   IKET=IVEC
        ! END IF

C Loop over types and symmetry block of VEC1 vector:
      DO 400 ICASE1=1,13
        DO 401 ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) GOTO 401
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          NIN1=NINDEP(ISYM1,ICASE1)
          NVEC1=NIS1*NAS1
          IF(NVEC1.EQ.0) GOTO 401
C Form VEC1 from the BRA vector, transformed to covariant form.
          !! Prepare vectors
          IMLTOP=1
          !! IBRA)
          Call PrepVec1(IMLTOP,NAS1,NIS1,ICASE1,ISYM1,NWEC1,
     *                  LVEC1,LWEC1,LVEC1S,LWEC1S,IVEC,ILOOP.EQ.1)
          If (ICASE1.LE.11) Then
            idSDer = idSDMat(iSym1,iCase1)
            CALL DDAFILE(LuSTD,2,Work(ipSDER1),nAS1*nAS1,idSDer)
          End If
C
          DO 300 ICASE2=ICASE1+1,13
            IF(IFCOUP(ICASE2,ICASE1).EQ.0) GOTO 300
C           if (icase1.ne.10.or.icase2.ne.12) cycle
C           if (icase1.ne. 8.or.icase2.ne.10) cycle
C           if (icase1.ne. 6.or.icase2.ne.12) cycle
            DO 200 ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) GOTO 200
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NIN2=NINDEP(ISYM2,ICASE2)
              NVEC2=NIS2*NAS2
              IF(NVEC2.EQ.0) GOTO 200
C
              IMLTOP=1
              !! IKET)
              Call PrepVec1(IMLTOP,NAS2,NIS2,ICASE2,ISYM2,NWEC2,
     *                      LVEC2,LWEC2,LVEC2S,LWEC2S,IVEC,ILOOP.EQ.2)
C
              !! S1*C1 derivative
              If (iCase1.LE.11) Call C1S1DER(Work(ipSDER1))
              !! C2 derivative
              If (iCase2.LE.11) Then
                idSDer = idSDMat(iSym2,iCase2)
                CALL DDAFILE(LuSTD,2,Work(ipSDER2),nAS2*nAS2,idSDer)
                Call C2DER(Work(ipSDER2))
                idSDer = idSDMat(iSym2,iCase2)
                CALL DDAFILE(LuSTD,1,Work(ipSDER2),nAS2*nAS2,idSDer)
              End If
C
              Call PrepVec2(NAS2,NIS2,ICASE2,NWEC2,
     *                      LVEC2,LWEC2,LVEC2S,LWEC2S)
C             CALL RHS_FREE(NAS2,NIS2,LVEC2)
 200        CONTINUE
 300      CONTINUE
          Call PrepVec2(NAS1,NIS1,ICASE1,NWEC1,
     *                  LVEC1,LWEC1,LVEC1S,LWEC1S)
          If (iCase1.LE.11) Then
            idSDer = idSDMat(iSym1,iCase1)
            CALL DDAFILE(LuSTD,1,Work(ipSDER1),nAS1*nAS1,idSDer)
          End If
C         CALL RHS_FREE(NAS1,NIS1,LVEC1)
C         IF(NWEC1.GT.0)
C    &         CALL GETMEM('WEC1','FREE','REAL',LWEC1,NWEC1)
 401    CONTINUE
 400  CONTINUE

 1000 CONTINUE
C
C
C
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSGM=CPUSGM+(CPU1-CPU0)
      TIOSGM=TIOSGM+(TIO1-TIO0)
C
      Call GETMEM('SDER1','FREE','REAL',ipSDER1,MaxLen)
      Call GETMEM('SDER2','FREE','REAL',ipSDER2,MaxLen)
C
      Call mma_deallocate(FIT_Full)
      Call mma_deallocate(FTI_Full)
      Call mma_deallocate(FIA_Full)
      Call mma_deallocate(FAI_Full)
      Call mma_deallocate(FTA_Full)
      Call mma_deallocate(FAT_Full)
      Do iSym = 1, nSym
         FIT(iSym)%A => Null()
         FTI(iSym)%A => Null()
         FIA(iSym)%A => Null()
         FAI(iSym)%A => Null()
         FTA(iSym)%A => Null()
         FAT(iSym)%A => Null()
      End Do

C Transform contrav C  to eigenbasis of H0(diag):
      CALL PTRTOSR(1,IVEC,IVEC)
      IF(IVEC.NE.JVEC) CALL PTRTOSR(1,JVEC,JVEC)

C 99  CONTINUE
      RETURN

      Contains
C
C-----------------------------------------------------------------------
C
      Subroutine PrepVec1(IMLTOP_,nAS_,nIS_,iCase_,iSym_,nWec_,
     *                    lVec_,lWec_,lVecS_,lWecS_,iVec_,ADDLAM)
C
      Implicit Real*8 (A-H,O-Z)
C
      Logical ADDLAM
C
      !! Prepare VecS and WecS
      CALL RHS_ALLO(nAS_,nIS_,lVec_)
      CALL RHS_READ(nAS_,nIS_,lVec_,iCase_,iSym_,iVec_)
      If (IVEC.NE.JVEC.and.ADDLAM) Then
        !! T = T + \lambda
        If (SCAL.ne.1.0D+00) Call DScal_(nAS_*nIS_,SCAL,Work(lVec_),1)
        CALL RHS_ALLO(nAS_,nIS_,LTMP)
        CALL RHS_READ(nAS_,nIS_,LTMP,iCase_,iSym_,JVEC)
        Call DaXpY_(nAS_*nIS_,1.0D+00,Work(LTMP),1,Work(lVec_),1)
        CALL RHS_FREE(nAS_,nIS_,LTMP)
      End If
      nWec_ = 0
      lWec_ = 1
      FACT  = 1.0D00/(DBLE(MAX(1,NACTEL)))
      IF(iCase_.EQ.1) nWec_ =NASH(iSym_)*NISH(iSym_)
      IF(iCase_.EQ.4) nWec_ =NASH(iSym_)*NSSH(iSym_)
      IF(iCase_.EQ.5.AND.iSym_.EQ.1) nWec_ = nIS_
      !! Consider parallelization later...
      IF(nWec_.GT.0) THEN
        CALL GETMEM('WEC1','ALLO','REAL',lWec_,nWec_)
        CALL DCOPY_(nWec_,[0.0D0],0,WORK(lWec_),1)
        IF(iCase_.EQ.1) THEN
          CALL SPEC1A(IMLTOP_,FACT,iSym_,WORK(lVec_),WORK(lWec_))
        ELSE IF(iCase_.EQ.4) THEN
          CALL SPEC1C(IMLTOP_,FACT,iSym_,WORK(lVec_),WORK(lWec_))
        ELSE IF(iCase_.EQ.5.AND.iSym_.EQ.1) THEN
          CALL SPEC1D(IMLTOP_,FACT,WORK(lVec_),WORK(lWec_))
        END IF
      END IF
C
C
C
      !! is it OK?
      If (iCase_.GT.11) Then
        lVecS_ = lVec_
        lWecS_ = lWec_
        Return
      End If
C
C
C
      !! Prepare VecS and WecS
      CALL RHS_ALLO(nAS_,nIS_,lVecS_)
      CALL RHS_SCAL(nAS_,nIS_,lVecS_,0.0D+00)
      CALL RHS_STRANS (nAS_,nIS_,1.0D+00,lVec_,lVecS_,iCase_,iSym_)
C
      IF(nWec_.GT.0) THEN
        CALL GETMEM('WEC1S','ALLO','REAL',lWecS_,nWec_)
        CALL DCOPY_(nWec_,[0.0D0],0,WORK(lWecS_),1)
        IF(iCase_.EQ.1) THEN
          CALL SPEC1A(IMLTOP_,FACT,iSym_,WORK(lVecS_),WORK(lWecS_))
        ELSE IF(iCase_.EQ.4) THEN
          CALL SPEC1C(IMLTOP_,FACT,iSym_,WORK(lVecS_),WORK(lWecS_))
        ELSE IF(iCase_.EQ.5.AND.iSym_.EQ.1) THEN
          CALL SPEC1D(IMLTOP_,FACT,WORK(lVecS_),WORK(lWecS_))
        END IF
      END IF
C
      End Subroutine PrepVec1
C
C-----------------------------------------------------------------------
C
      Subroutine PrepVec2(nAS_,nIS_,iCase_,nWec_,
     *                    lVec_,lWec_,lVecS_,lWecS_)
C
      Implicit Real*8 (A-H,O-Z)
C
      CALL RHS_FREE(nAS_,nIS_,lVec_)
C     If (nWec_.GT.0) CALL RHS_FREE(nAS_n,nIS_,lWec_)
      If (nWec_.GT.0) CALL GETMEM('WEC1','FREE','REAL',lWec_,nWec_)
C
      If (iCase_.GT.11) Return
C
      CALL RHS_FREE(nAS_,nIS_,lVecS_)
C     If (nWec_.GT.0) CALL RHS_FREE(nAS_n,nIS_,lWecS_)
      If (nWec_.GT.0) CALL GETMEM('WEC1S','FREE','REAL',lWecS_,nWec_)
C
      End Subroutine PrepVec2
C
C-----------------------------------------------------------------------
C
      Subroutine C1S1DER(SDER)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension SDER(*)
C
C     (T2Ct2*f)py * (T1Ct1)pz * dS1yz/da
C     -1/2 (T2Ct2*f*S1*C1*Ct1)pt * (T1Ct1)pu * dS1tu/da
C
      !! initialize
      CALL GETMEM('TMP2','ALLO','REAL',LTMP2,NVEC1)
      CALL DCOPY_(NVEC1,[0.0D0],0,WORK(LTMP2),1)
      CALL GETMEM('TMP1','ALLO','REAL',LTMP1,MAX(1,NWEC1))
      IF(NWEC1.GT.0) THEN
        CALL DCOPY_(NWEC1,[0.0D0],0,WORK(LTMP1),1)
      END IF
C
      !! 1. T2*Ct2*f
      IMLTOP=0
      CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &         WORK(LTMP1),LTMP2,LVEC2,iWORK(LLISTS))
C
      IF(NWEC1.GT.0) THEN
        FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
        IF (ICASE1.EQ.1) THEN
          CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LTMP2),
     &              WORK(LTMP1))
        ELSE IF(ICASE1.EQ.4) THEN
          CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LTMP2),
     &              WORK(LTMP1))
        ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
          CALL SPEC1D(IMLTOP,FACT,WORK(LTMP2),WORK(LTMP1))
        END IF
      END IF
      CALL GETMEM('TMP1','FREE','REAL',LTMP1,MAX(1,NWEC1))
C
      !! Finalize the derivative of S1
      !! 2S. (T2Ct2*f) * T1Ct1
      Call DGEMM_('N','T',NAS1,NAS1,NIS1,
     *            2.0D+00,WORK(LVEC1),NAS1,WORK(LTMP2),NAS1,
     *            1.0D+00,SDER,NAS1)
C
      !! Next, the derivative of C1
      !! 2C. (T2Ct2*f) * S1*C1 (MO -> IC)
      CALL RHS_ALLO(NIN1,NIS1,LTMP1)
      ITYPE=1
      CALL RHS_SR2C (ITYPE,1,NAS1,NIS1,NIN1,LTMP1,LTMP2,ICASE1,ISYM1)
      !! 3C. (T2Ct2*f) * S1*C1 * Ct1 (IC -> MO)
      ITYPE=0
      CALL RHS_SR2C (ITYPE,0,NAS1,NIS1,NIN1,LTMP1,LTMP2,ICASE1,ISYM1)
      CALL RHS_FREE(NIN1,NIS1,LTMP1)
C
      !! 4C. (T1Ct1*f) * (T2Ct2St2*f*C1*Ct1)
      Call DGEMM_('N','T',NAS1,NAS1,NIS1,
     *           -1.0D+00,WORK(LVEC1),NAS1,WORK(LTMP2),NAS1,
     *            1.0D+00,SDER,NAS1)
C
      CALL GETMEM('TMP2','FREE','REAL',LTMP2,NVEC1)
C
      End Subroutine C1S1DER
C
C-----------------------------------------------------------------------
C
      Subroutine C2DER(SDER)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension SDER(*)
C
C     -1/2 (T2Ct2)pu * dS2tu/da * (T1Ct1St1*f*C2*Ct2)pt
C
      !! initialize
      CALL RHS_ALLO(NAS2,NIS2,LTMP)
      CALL RHS_SCAL(NAS2,NIS2,LTMP,0.0D+00)
C
      !! 1. T1*Ct1*St1*f
      IMLTOP=1
      CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &         WORK(LWEC1S),LVEC1S,LTMP,iWORK(LLISTS))
C
      !! 2. (T1Ct1St1*f) * C2 (MO -> IC)
      CALL RHS_ALLO(NIN2,NIS2,LTMP2)
      ITYPE=0
      CALL RHS_SR2C (ITYPE,1,NAS2,NIS2,NIN2,LTMP2,LTMP,ICASE2,ISYM2)
      !! 3. (T1Ct1St1*f) * C2 * Ct2 (IC -> MO)
      CALL RHS_SR2C (ITYPE,0,NAS2,NIS2,NIN2,LTMP2,LTMP,ICASE2,ISYM2)
      CALL RHS_FREE(NIN2,NIS2,LTMP2)
C
      !! 4. (T2Ct2*f) * (T1Ct1St1*f*C2*Ct2)
      Call DGEMM_('N','T',NAS2,NAS2,NIS2,
     *           -1.0D+00,WORK(LVEC2),NAS2,WORK(LTMP),NAS2,
     *            1.0D+00,SDER,NAS2)
C
      CALL RHS_FREE(NAS2,NIS2,LTMP)
C
      End Subroutine C2DER
C
      End subroutine sigder
