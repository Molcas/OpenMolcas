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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VSCT(*)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Read coefficient vector from LUSOLV (C repres).
#ifdef _DEBUG_
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
      IDS=IWORK(LIDSCT-1+ISCT+MXSCT*(ISYM-1+8*
     &                         (ICASE-1+MXCASE*(IVEC-1))))
      NCOL=MIN(NIS-MDVEC*(ISCT-1),MDVEC)
      NSCT=NAS*NCOL
      CALL DDAFILE(LUSOLV,2,VSCT,NSCT,IDS)
#ifdef _DEBUG_
        WRITE(6,*)' First few elements:'
        WRITE(6,'(1x,5f15.6)')(VSCT(I),I=1,MIN(NSCT,10))
#endif
      RETURN
      END

      SUBROUTINE RDBLKC(ISYM,ICASE,IVEC,VEC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VEC(*)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Read coefficient vector from LUSOLV (C repres).
#ifdef _DEBUG_
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
      IDV=IWORK(LIDSCT+MXSCT*(ISYM-1+8*
     &                         (ICASE-1+MXCASE*(IVEC-1))))
      LVEC=1
      DO 10 IISTA=1,NIS,MDVEC
        NCOL=MIN(NIS+1-IISTA,MDVEC)
        NBLK=NAS*NCOL
        CALL DDAFILE(LUSOLV,2,VEC(LVEC),NBLK,IDV)
        LVEC=LVEC+NBLK
  10  CONTINUE
#ifdef _DEBUG_
        WRITE(6,*)' First few elements:'
        WRITE(6,'(1x,5f15.6)')(VEC(I),I=1,MIN(NCOEF,10))
#endif
      RETURN
      END

************************************************************************
* New VECTOR UTILITIES, written by Steven Vancoillie, May 2011
* A set of subroutines that can transform RHS arrays using the parallel
* aware subroutines
************************************************************************

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PSCAVEC (FACT,IVEC,JVEC)
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

      REAL*8 FACT
      INTEGER IVEC,JVEC

      REAL*8 CPU0,CPU1,CPU,TIO0,TIO1,TIO
      INTEGER ICASE,ISYM,NIN,NIS
      INTEGER lg_V

C Scale vector nr IVEC with scale factor FACT and put the result in
C vector nr JVEC: |JVEC> <- FACT * |IVEC>
      CALL TIMING(CPU0,CPU,TIO0,TIO)

      IF(FACT.EQ.1.0D00.AND.IVEC.EQ.JVEC) RETURN
      DO ICASE=1,NCASES
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NIN*NIS.NE.0) THEN
            CALL RHS_ALLO (NIN,NIS,lg_V)
            CALL RHS_READ (NIN,NIS,lg_V,ICASE,ISYM,IVEC)
            CALL RHS_SCAL (NIN,NIS,lg_V,FACT)
            CALL RHS_SAVE (NIN,NIS,lg_V,ICASE,ISYM,JVEC)
            CALL RHS_FREE (NIN,NIS,lg_V)
          END IF
        END DO
      END DO

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSCA=CPUSCA+(CPU1-CPU0)
      TIOSCA=TIOSCA+(TIO1-TIO0)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE POVLVEC (IVEC,JVEC,OVLAPS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

      DIMENSION OVLAPS(0:8,0:MXCASE)

C Compute overlaps of vectors nr IVEC and JVEC in SR format!, for each
C individual case and symmetry block, in OVLAPS(ISYM,ICASE), summed over
C symmetry in OVLAPS(0,ICASE), summed over case in OVLAPS(ISYM,0), total
C sum in OVLAPS(0,0).
      CALL TIMING(CPU0,CPU,TIO0,TIO)

      OVLTOT=0.0D0
      DO ISYM=1,NSYM
        OVLAPS(ISYM,0)=0.0D0
      END DO
      DO ICASE=1,NCASES
        OVLSUM=0.0D0
        DO ISYM=1,NSYM
          OVL=0.0D0
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
            CALL RHS_FREE (NIN,NIS,lg_V1)
            IF (IVEC.NE.JVEC) THEN
              CALL RHS_FREE (NIN,NIS,lg_V2)
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

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PLCVEC (ALPHA,BETA,IVEC,JVEC)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

C |JVEC> := BETA*|JVEC> + ALPHA*|IVEC>, IVEC and JVEC in SR format!

      IF(BETA.EQ.1.0D0.AND.ALPHA.EQ.0.0D0) RETURN

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      DO 100 ICASE=1,NCASES
        DO 101 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIN*NIS.EQ.0) GOTO 101

          CALL RHS_ALLO (NIN,NIS,lg_V2)
          IF(BETA.NE.0.0D0.AND.ALPHA.NE.0.0D0) THEN
            CALL RHS_ALLO (NIN,NIS,lg_V1)
          END IF

          IF(BETA.EQ.0.0D0.AND.ALPHA.EQ.0.0D0) THEN
              CALL RHS_SCAL (NIN,NIS,lg_V2,0.0D0)
          ELSE IF(BETA.NE.0.0D0) THEN
            CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
            IF (ALPHA.NE.0.0D0) THEN
              CALL RHS_SCAL (NIN,NIS,lg_V2,BETA)
              CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
              CALL RHS_DAXPY (NIN,NIS,ALPHA,lg_V1,lg_V2)
            ELSE
              CALL RHS_SCAL (NIN,NIS,lg_V2,BETA)
            END IF
          ELSE IF(ALPHA.NE.0.0D0) THEN
            CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,IVEC)
            CALL RHS_SCAL (NIN,NIS,lg_V2,ALPHA)
          END IF

          CALL RHS_SAVE (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)

          CALL RHS_FREE (NIN,NIS,lg_V2)
          IF(ALPHA.NE.0.0D0) THEN
            CALL RHS_FREE (NIN,NIS,lg_V1)
          END IF
 101    CONTINUE
 100  CONTINUE

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPULCS=CPULCS+(CPU1-CPU0)
      TIOLCS=TIOLCS+(TIO1-TIO0)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PTRTOC (ITYPE,IVEC,JVEC)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

C Transform RHS vectors from SR format to C format.
C ITYPE=0 uses only T matrix, ITYPE=1 uses S*T matrix

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      DO 200 ICASE=1,NCASES
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.12) GOTO 200
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.13) GOTO 200
        DO 100 ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          IF(NAS.EQ.0) GOTO 100
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) GOTO 100
          CALL RHS_ALLO (NAS,NIS,lg_V2)
          IF(ICASE.NE.12 .AND. ICASE.NE.13) THEN
            IF(NIN.GT.0) THEN
              CALL RHS_ALLO (NIN,NIS,lg_V1)
              CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
              CALL RHS_SR2C (ITYPE,0,NAS,NIS,NIN,lg_V1,lg_V2,ICASE,ISYM)
              CALL RHS_FREE (NIN,NIS,lg_V1)
            ELSE
              CALL RHS_SCAL (NAS,NIS,lg_V2,0.0D0)
            END IF
          ELSE
            CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,IVEC)
          END IF

          CALL RHS_SAVE (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)

          CALL RHS_FREE (NAS,NIS,lg_V2)

 100    CONTINUE
 200  CONTINUE

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUVEC=CPUVEC+(CPU1-CPU0)
      TIOVEC=TIOVEC+(TIO1-TIO0)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE PTRTOSR (ITYPE,IVEC,JVEC)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

C Transform RHS vectors from SR format to C format.
C ITYPE=0 uses only T matrix, ITYPE=1 uses S*T matrix

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      DO 200 ICASE=1,NCASES
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.12) GOTO 200
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.13) GOTO 200
        DO 100 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 100
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) GOTO 100
          CALL RHS_ALLO (NIN,NIS,lg_V1)
          IF(ICASE.NE.12 .AND. ICASE.NE.13) THEN
            IF(NAS.GT.0) THEN
              CALL RHS_ALLO (NAS,NIS,lg_V2)
              CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,IVEC)
              CALL RHS_SR2C (ITYPE,1,NAS,NIS,NIN,lg_V1,lg_V2,ICASE,ISYM)
              CALL RHS_FREE (NAS,NIS,lg_V2)
            ELSE
              CALL RHS_SCAL (NIN,NIS,lg_V1,0.0D0)
            END IF
          ELSE
            CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
          END IF

          CALL RHS_SAVE (NIN,NIS,lg_V1,ICASE,ISYM,JVEC)

          CALL RHS_FREE (NIN,NIS,lg_V1)

 100    CONTINUE
 200  CONTINUE

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUVEC=CPUVEC+(CPU1-CPU0)
      TIOVEC=TIOVEC+(TIO1-TIO0)

      RETURN
      END
