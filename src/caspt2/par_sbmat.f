************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*
* WRAPPERS FOR PARALLEL S AND B MATRIX ROUTINES
*
      SUBROUTINE PSBMAT_GETMEM(cNAME,lg_M,nSize)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_ga, only: GA_arrays, Allocate_GA_Array
      IMPLICIT None
CSVC2010: create square global array S/B for symmetry iSYM
C with integer handle lg_M or if replicate or serial, create
C tridiagonal local array at Work(lg_M)
#include "caspt2.fh"
#include "pt2_guga.fh"
      Integer lg_M, nSize
      CHARACTER(len=*) cNAME

      Integer nTri

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_CREATE_STRIPED ('H',nSize,nSize,cNAME,LG_M)
        CALL GA_ZERO (LG_M)
      ELSE
#endif
        nTri=(nSize*(nSize+1))/2
        lg_M=Allocate_GA_Array(nTri,cName)
        GA_Arrays(lg_M)%A(:)=0.0D0
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE PSBMAT_GETMEM

      SUBROUTINE PSBMAT_FREEMEM(lg_M)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_ga, only: Deallocate_GA_Array
      IMPLICIT NONE
CSVC2010: destroy square global array S/B for symmetry iSYM
C with integer handle lg_M or if replicate or serial, free the
C tridiagonal local array at Work(lg_M)
#include "caspt2.fh"
#include "pt2_guga.fh"
      Integer lg_M

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL bStat
#endif


#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        bStat = GA_Destroy(lg_M)
      ELSE
#endif
        Call Deallocate_GA_Array(lg_M)
#ifdef _MOLCAS_MPP_
      END IF
#endif
      END SUBROUTINE PSBMAT_FREEMEM

      SUBROUTINE PSBMAT_WRITE(cNAME,iCase,iSym,lg_M,nSize)
CSVC20100902: write the global array lg_M to disk using DRA interface,
C or if replicate or serial, write WORK(lg_M) to LUSBT
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use caspt2_global, only: LUH0T
#endif
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT, IDBMAT, IDTMAT, IDSTMAT
      use fake_ga, only: GA_arrays
      IMPLICIT None
#include "caspt2.fh"
#include "pt2_guga.fh"
      Integer iCase, iSym, lg_M, nSize
      CHARACTER(LEN=*) cNAME


#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      Integer LU, myRank, ISTA,IEND,JSTA,JEND, mpt_M, LDM
#endif
      Integer IDISK, nBlock


      IF (CNAME.EQ.'S') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(1)
#endif
        IDISK=IDSMAT(iSym,iCase)
        nBlock=(nSize*(nSize+1))/2
      ELSE IF (CNAME.EQ.'B') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(2)
#endif
        IDISK=IDBMAT(iSym,iCase)
        nBlock=(nSize*(nSize+1))/2
      ELSE IF (CNAME.EQ.'T') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(3)
#endif
        IDISK=IDTMAT(iSym,iCase)
        nBlock=nSize
      ELSE IF (CNAME.EQ.'M') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(4)
#endif
        IDISK=IDSTMAT(iSym,iCase)
        nBlock=nSize
      END IF

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_M,myRank,ISTA,IEND,JSTA,JEND)
        IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
          CALL GA_Access (lg_M,ISTA,IEND,JSTA,JEND,mpt_M,LDM)
          NBLOCK=LDM*(JEND-JSTA+1)
          CALL DDAFILE(LU,1,DBL_MB(mpt_M),NBLOCK,IDISK)
          CALL GA_Release (lg_M,ISTA,IEND,JSTA,JEND)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        CALL DDAFILE(LUSBT,1,GA_Arrays(lg_M)%A(:),nBlock,IDISK)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE PSBMAT_WRITE

      SUBROUTINE PSBMAT_READ(cNAME,iCase,iSym,lg_M,nSize)
CSVC20100902: read the disk array stored as cName+iSym using DRA
C interface into global array lg_M, or if replicate or serial, read from
C LUSBT into WORK(lg_M)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use caspt2_global, only: LUH0T
#endif
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT, IDBMAT, IDTMAT, IDSTMAT
      use fake_ga, only: GA_arrays
      IMPLICIT None
#include "caspt2.fh"
#include "pt2_guga.fh"
      INTEGER iCASE,iSym,lg_M,nSize
      CHARACTER(LEN=*) cNAME

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      INTEGER LU, myRank, iSTA, IEND, JSTA, JEND, mpt_M, LDM
#endif
      INTEGER IDISK, nBlock


      IF (CNAME.EQ.'S') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(1)
#endif
        IDISK=IDSMAT(iSym,iCase)
        nBlock=(nSize*(nSize+1))/2
      ELSE IF (CNAME.EQ.'B') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(2)
#endif
        IDISK=IDBMAT(iSym,iCase)
        nBlock=(nSize*(nSize+1))/2
      ELSE IF (CNAME.EQ.'T') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(3)
#endif
        IDISK=IDTMAT(iSym,iCase)
        nBlock=nSize
      ELSE IF (CNAME.EQ.'M') THEN
#ifdef _MOLCAS_MPP_
        LU=LUH0T(4)
#endif
        IDISK=IDSTMAT(iSym,iCase)
        nBlock=nSize
      END IF

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_M,myRank,ISTA,IEND,JSTA,JEND)
        IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
          CALL GA_Access (lg_M,ISTA,IEND,JSTA,JEND,mpt_M,LDM)
          NBLOCK=LDM*(JEND-JSTA+1)
          CALL DDAFILE(LU,2,DBL_MB(mpt_M),NBLOCK,IDISK)
          CALL GA_Release_Update (lg_M,ISTA,IEND,JSTA,JEND)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        CALL DDAFILE(LUSBT,2,GA_Arrays(lg_M)%A(:),nBlock,IDISK)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE PSBMAT_READ

      REAL*8 FUNCTION PSBMAT_FPRINT(lg_M,NM)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_ga, only: GA_arrays
      IMPLICIT NONE
      INTEGER lg_M, NM

      INTEGER nTri
      REAL*8, External::DNRM2_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        PSBMAT_FPRINT=SQRT(GA_DDOT(lg_M,lg_M))
      ELSE
#endif
        nTri=(NM*(NM+1))/2
        PSBMAT_FPRINT=DNRM2_(nTri,GA_Arrays(lg_M)%A(:),1)
#ifdef _MOLCAS_MPP_
      END IF
#endif
      END FUNCTION PSBMAT_FPRINT
