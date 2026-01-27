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
* Copyright (C) 2008, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2008  PER-AKE MALMQUIST                    *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE RDSCTC(ISCT,ISYM,ICASE,IVEC,VSCT)
      use definitions, only: iwp, wp
      use caspt2_global, only: LUSOLV, IDSCT
      use EQSOLV, only: MxSCT, ModVec
      use caspt2_module, only: NASUP, NISUP, MxCASE
#ifdef _DEBUGPRINT_
      use caspt2_module, only: cases
#endif
      IMPLICIT NONE
      real(kind=wp), intent(out) :: VSCT(*)
      integer(kind=iwp), intent(In):: ISCT, iSYM, ICASE, IVEC

      integer(kind=iwp) iDS, MDVEC, NAS, NIS, NCOEF, NCOL, NSCT
#ifdef _DEBUGPRINT_
      integer(kind=iwp) I
#endif


C Read coefficient vector from LUSOLV (C repres).
#ifdef _DEBUGPRINT_
        WRITE(6,*)' RDSCTC (Normal repres.)'
        WRITE(6,'(a,i2,a,a,a,i2,a,i2)')' Vector nr.',IVEC,
     &          '  Case ',CASES(ICASE),' Symm ',ISYM,
     &          ' Section ',ISCT
#endif
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      NCOEF=NAS*NIS
      IF(NCOEF.EQ.0) RETURN
      MDVEC=MODVEC(ISYM,ICASE)
*      IDS=IDSCT(ISCT,ISYM,ICASE,IVEC)
      IDS=IDSCT(ISCT+MXSCT*(ISYM-1+8*
     &                         (ICASE-1+MXCASE*(IVEC-1))))
      NCOL=MIN(NIS-MDVEC*(ISCT-1),MDVEC)
      NSCT=NAS*NCOL
      CALL DDAFILE(LUSOLV,2,VSCT,NSCT,IDS)
#ifdef _DEBUGPRINT_
        WRITE(6,*)' First few elements:'
        WRITE(6,'(1x,5f15.6)')(VSCT(I),I=1,MIN(NSCT,10))
#endif
      END SUBROUTINE RDSCTC

      SUBROUTINE RDBLKC(ISYM,ICASE,IVEC,VEC)
      use definitions, only: iwp, wp
      use caspt2_global, only: LUSOLV, IDSCT
      use EQSOLV, only: MxSct, ModVec
      use caspt2_module, only: NASUP, NISUP, MxCase
#ifdef _DEBUGPRINT_
      use caspt2_module, only: CASES
#endif
      IMPLICIT None
      real(kind=wp), intent(out):: VEC(*)
      integer(kind=iwp), intent(in):: iSym,iCase,iVec

      integer(kind=iwp) NAS, NIS, NCOEF, MDVEC, IDV, LVEC, IISTA,
     &                  NCOL, NBLK
#ifdef _DEBUGPRINT_
      integer(kind=iwp) I
#endif

C Read coefficient vector from LUSOLV (C repres).
#ifdef _DEBUGPRINT_
        WRITE(6,*)' RDBLKC (Normal repres.)'
        WRITE(6,'(a,i2,a,a,a,i2)')' Vector nr.',IVEC,
     &          '  Case ',CASES(ICASE),' Symm ',ISYM
#endif
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      NCOEF=NAS*NIS
      IF(NCOEF.EQ.0) RETURN
      MDVEC=MODVEC(ISYM,ICASE)
*      IDV=IDSCT(1,ISYM,ICASE,IVEC)
      IDV=IDSCT(1+MXSCT*(ISYM-1+8*
     &                         (ICASE-1+MXCASE*(IVEC-1))))
      LVEC=1
      DO IISTA=1,NIS,MDVEC
        NCOL=MIN(NIS+1-IISTA,MDVEC)
        NBLK=NAS*NCOL
        CALL DDAFILE(LUSOLV,2,VEC(LVEC),NBLK,IDV)
        LVEC=LVEC+NBLK
      End Do
#ifdef _DEBUGPRINT_
        WRITE(6,*)' First few elements:'
        WRITE(6,'(1x,5f15.6)')(VEC(I),I=1,MIN(NCOEF,10))
#endif
      END SUBROUTINE RDBLKC

************************************************************************
* New VECTOR UTILITIES, written by Steven Vancoillie, May 2011
* A set of subroutines that can transform RHS arrays using the parallel
* aware subroutines
************************************************************************

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PSCAVEC (FACT,IVEC,JVEC)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_global, ONLY: iPrGlb
      USE PrintLevel, ONLY: usual
      use caspt2_module, only: CPUSCA, nCases, nSym, TIOSCA, nInDep,
     &                         niSup

      IMPLICIT NONE

      real(kind=wp), intent(in) ::  FACT
      integer(kind=iwp), intent(in) :: IVEC,JVEC

      real(kind=wp) :: CPU0,CPU1,CPU,TIO0,TIO1,TIO
      integer(kind=iwp) :: ICASE,ISYM,NIN,NIS
      integer(kind=iwp) :: lg_V

      REAL(kind=wp) ::  SIGMA2
      real(kind=wp), EXTERNAL :: RHS_DDOT

C Scale vector nr IVEC with scale factor FACT and put the result in
C vector nr JVEC: |JVEC> <- FACT * |IVEC>
      CALL TIMING(CPU0,CPU,TIO0,TIO)

      IF(FACT.EQ.One.AND.IVEC.EQ.JVEC) RETURN
      SIGMA2=Zero
      DO ICASE=1,NCASES
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NIN*NIS.NE.0) THEN
            CALL RHS_ALLO (NIN,NIS,lg_V)
            CALL RHS_READ (NIN,NIS,lg_V,ICASE,ISYM,IVEC)
            CALL RHS_SCAL (NIN,NIS,lg_V,FACT)
            IF(FACT.EQ.-One)SIGMA2=SIGMA2
     &                               +RHS_DDOT(NIN,NIS,lg_V,lg_V)
            CALL RHS_SAVE (NIN,NIS,lg_V,ICASE,ISYM,JVEC)
            CALL RHS_FREE (lg_V)
          END IF
        END DO
      END DO
      IF ((IPRGLB.GE.USUAL).AND.(FACT.EQ.-One)) THEN ! it is at ITER=0
         WRITE(6,*)
         WRITE(6,'(1x,a,f18.10)') 'Variance of |WF0>: ',SIGMA2
      END IF

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSCA=CPUSCA+(CPU1-CPU0)
      TIOSCA=TIOSCA+(TIO1-TIO0)

      END SUBROUTINE PSCAVEC

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE POVLVEC (IVEC,JVEC,OVLAPS)
      use definitions, only: wp, iwp
      use constants, only: Zero
      use caspt2_module, only: NINDEP, NISUP, MxCase, nCASES, CPUOVL,
     &                         TIOOVL, nSym
      IMPLICIT None

      real(kind=wp), intent(out) :: OVLAPS(0:8,0:MXCASE)
      integer(kind=iwp), intent(in) :: iVec, jVec

      real(kind=wp) CPU, CPU0, CPU1
      real(kind=wp) TIO, TIO0, TIO1
      real(kind=wp) OVLTOT, OVL, OVLSUM
      integer(kind=iwp) iCase, iSym, lg_v1, lg_v2, NIN, NIS
      real(kind=wp), External :: RHS_DDOT

C Compute overlaps of vectors nr IVEC and JVEC in SR format!, for each
C individual case and symmetry block, in OVLAPS(ISYM,ICASE), summed over
C symmetry in OVLAPS(0,ICASE), summed over case in OVLAPS(ISYM,0), total
C sum in OVLAPS(0,0).
      CALL TIMING(CPU0,CPU,TIO0,TIO)

      OVLTOT=Zero
      OVLAPS(:,0)=Zero

      DO ICASE=1,NCASES
        OVLSUM=Zero
        DO ISYM=1,NSYM
          OVL=Zero
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NIN*NIS.NE.0) THEN
            CALL RHS_ALLO (NIN,NIS,lg_V1)
            CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
            IF (IVEC.NE.JVEC) THEN
              CALL RHS_ALLO (NIN,NIS,lg_V2)
              CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
            ELSE
              lg_V2=lg_V1
            END IF
            OVL=RHS_DDOT(NIN,NIS,lg_V1,lg_V2)
            CALL RHS_FREE (lg_V1)
            IF (IVEC.NE.JVEC) THEN
              CALL RHS_FREE (lg_V2)
            END IF
          END IF
          OVLAPS(ISYM,ICASE)=OVL
          OVLAPS(ISYM,0)=OVLAPS(ISYM,0)+OVL
          OVLSUM=OVLSUM+OVL
        END DO
        OVLAPS(0,ICASE)=OVLSUM
        OVLTOT=OVLTOT+OVLSUM
      END DO
      OVLAPS(0,0)=OVLTOT

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUOVL=CPUOVL+(CPU1-CPU0)
      TIOOVL=TIOOVL+(TIO1-TIO0)

      END SUBROUTINE POVLVEC

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PLCVEC (ALPHA,BETA,IVEC,JVEC)
      use constants, only: Zero, One
      use definitions, only: iwp, wp
      use caspt2_module, only: nCases, nSym, nInDep, NISUP, CPULCS,
     &                         TIOLCS
      IMPLICIT None

      real(kind=wp), intent(in):: Alpha, Beta
      integer(kind=iwp), intent(in):: iVec, jVec

      real(kind=wp) CPU, CPU0, CPU1
      real(kind=wp) TIO, TIO0, TIO1
      integer(kind=iwp) iCase, iSym, lg_v1, lg_v2
      integer(kind=iwp) NIN, NIS

C |JVEC> := BETA*|JVEC> + ALPHA*|IVEC>, IVEC and JVEC in SR format!

      IF(BETA==One.AND.ALPHA==Zero) RETURN

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      DO ICASE=1,NCASES
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIN*NIS.EQ.0) Cycle

          CALL RHS_ALLO (NIN,NIS,lg_V2)
          IF(BETA.NE.Zero.AND.ALPHA.NE.Zero) THEN
            CALL RHS_ALLO (NIN,NIS,lg_V1)
          END IF

          IF(BETA.EQ.Zero.AND.ALPHA.EQ.Zero) THEN
              CALL RHS_SCAL (NIN,NIS,lg_V2,Zero)
          ELSE IF(BETA.NE.Zero) THEN
            CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
            IF (ALPHA.NE.Zero) THEN
              CALL RHS_SCAL (NIN,NIS,lg_V2,BETA)
              CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
              CALL RHS_DAXPY (NIN,NIS,ALPHA,lg_V1,lg_V2)
            ELSE
              CALL RHS_SCAL (NIN,NIS,lg_V2,BETA)
            END IF
          ELSE IF(ALPHA.NE.Zero) THEN
            CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,IVEC)
            CALL RHS_SCAL (NIN,NIS,lg_V2,ALPHA)
          END IF

          CALL RHS_SAVE (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)

          CALL RHS_FREE (lg_V2)
          IF(BETA.NE.Zero.AND.ALPHA.NE.Zero) THEN
            CALL RHS_FREE (lg_V1)
          END IF
        End Do
      End Do

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPULCS=CPULCS+(CPU1-CPU0)
      TIOLCS=TIOLCS+(TIO1-TIO0)

      END SUBROUTINE PLCVEC

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PTRTOC (ITYPE,IVEC,JVEC)
      use definitions, only: iwp, wp
      use Constants, only: Zero
      use caspt2_module, only: nCases, nSym, nInDep, NISUP, CPUVEC,
     &                         TIOVEC, NASUP
      IMPLICIT None

      integer(kind=iwp), intent(In):: iType, iVec, jVec

      real(kind=wp) CPU, CPU0, CPU1
      real(kind=wp) TIO, TIO0, TIO1
      integer(kind=iwp) iCase, iSym, lg_v1, lg_v2, NAS, NIN, NIS

C Transform RHS vectors from SR format to C format.
C ITYPE=0 uses only T matrix, ITYPE=1 uses S*T matrix

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      DO ICASE=1,NCASES
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.12) Cycle
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.13) Cycle
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          IF(NAS.EQ.0) Cycle
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) Cycle
          CALL RHS_ALLO (NAS,NIS,lg_V2)
          IF(ICASE.NE.12 .AND. ICASE.NE.13) THEN
            IF(NIN.GT.0) THEN
              CALL RHS_ALLO (NIN,NIS,lg_V1)
              CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
              CALL RHS_SR2C (ITYPE,0,NAS,NIS,NIN,
     &                       lg_V1,lg_V2,ICASE,ISYM)
              CALL RHS_FREE (lg_V1)
            ELSE
              CALL RHS_SCAL (NAS,NIS,lg_V2,Zero)
            END IF
          ELSE
            CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,IVEC)
          END IF

          CALL RHS_SAVE (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)

          CALL RHS_FREE (lg_V2)

        End Do
      End Do

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUVEC=CPUVEC+(CPU1-CPU0)
      TIOVEC=TIOVEC+(TIO1-TIO0)

      END SUBROUTINE PTRTOC

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PTRTOSR (ITYPE,IVEC,JVEC)
      use definitions, only: iwp, wp
      use Constants, only: Zero
      use caspt2_module, only: nSym, nInDep, nASup, nISup, nCases,
     &                         CPUVec, TIOVec
      IMPLICIT None

      integer(kind=iwp), intent(In):: iTYPE, iVec, jVec

      integer(kind=iwp) :: iCase, iSym, nIn, nAs, nIs
      integer(kind=iwp) :: lg_v1, lg_v2
      real(kind=wp) CPU1, CPU0, TIO1, TIO0, CPU, TIO


C Transform RHS vectors from SR format to C format.
C ITYPE=0 uses only T matrix, ITYPE=1 uses S*T matrix

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      DO ICASE=1,NCASES
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.12) Cycle
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.13) Cycle
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) Cycle
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) Cycle
          CALL RHS_ALLO (NIN,NIS,lg_V1)
          IF(ICASE.NE.12 .AND. ICASE.NE.13) THEN
            IF(NAS.GT.0) THEN
              CALL RHS_ALLO (NAS,NIS,lg_V2)
              CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,IVEC)
              CALL RHS_SR2C (ITYPE,1,NAS,NIS,NIN,
     &                       lg_V1,lg_V2,ICASE,ISYM)
              CALL RHS_FREE (lg_V2)
            ELSE
              CALL RHS_SCAL (NIN,NIS,lg_V1,Zero)
            END IF
          ELSE
            CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
          END IF

          CALL RHS_SAVE (NIN,NIS,lg_V1,ICASE,ISYM,JVEC)

          CALL RHS_FREE (lg_V1)

        End Do
      End Do

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUVEC=CPUVEC+(CPU1-CPU0)
      TIOVEC=TIOVEC+(TIO1-TIO0)

      END SUBROUTINE PTRTOSR
