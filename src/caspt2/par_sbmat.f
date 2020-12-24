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
      IMPLICIT REAL*8 (A-H,O-Z)
CSVC2010: create square global array S/B for symmetry iSYM
C with integer handle lg_M or if replicate or serial, create
C tridiagonal local array at Work(lg_M)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"

#include "SysDef.fh"

      CHARACTER(len=*) cNAME
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif


#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_CREATE_STRIPED ('H',nSize,nSize,cNAME,LG_M)
        CALL GA_ZERO (LG_M)
      ELSE
        nTri=(nSize*(nSize+1))/2
        CALL GETMEM(cNAME,'ALLO','REAL',lg_M,nTri)
        CALL DCOPY_(nTri,[0.0D0],0,WORK(lg_M),1)
      END IF
#else
      nTri=(nSize*(nSize+1))/2
      CALL GETMEM(cNAME,'ALLO','REAL',lg_M,nTri)
      CALL DCOPY_(nTri,[0.0D0],0,WORK(lg_M),1)
#endif

      END

      SUBROUTINE PSBMAT_FREEMEM(cNAME,lg_M,nSize)
      IMPLICIT REAL*8 (A-H,O-Z)
CSVC2010: destroy square global array S/B for symmetry iSYM
C with integer handle lg_M or if replicate or serial, free the
C tridiagonal local array at Work(lg_M)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"

#include "SysDef.fh"

      CHARACTER(len=*) cNAME
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL bStat
#endif


#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        bStat = GA_Destroy(lg_M)
      ELSE
        nTri=(nSize*(nSize+1))/2
        CALL GETMEM(cNAME,'FREE','REAL',lg_M,nTri)
      END IF
#else
      nTri=(nSize*(nSize+1))/2
      CALL GETMEM(cNAME,'FREE','REAL',lg_M,nTri)
#endif
      END

      SUBROUTINE PSBMAT_WRITE(cNAME,iCase,iSym,lg_M,nSize)
CSVC20100902: write the global array lg_M to disk using DRA interface,
C or if replicate or serial, write WORK(lg_M) to LUSBT
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"

#include "SysDef.fh"
      CHARACTER cNAME

#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif


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
*       nTri=(nSize*(nSize+1))/2
        CALL DDAFILE(LUSBT,1,WORK(lg_M),nBlock,IDISK)
      END IF
#else
*     nTri=(nSize*(nSize+1))/2
      CALL DDAFILE(LUSBT,1,WORK(lg_M),nBlock,IDISK)
#endif

      END

      SUBROUTINE PSBMAT_READ(cNAME,iCase,iSym,lg_M,nSize)
CSVC20100902: read the disk array stored as cName+iSym using DRA
C interface into global array lg_M, or if replicate or serial, read from
C LUSBT into WORK(lg_M)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"

#include "SysDef.fh"
      CHARACTER cNAME

#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif


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
*       nTri=(nSize*(nSize+1))/2
        CALL DDAFILE(LUSBT,2,WORK(lg_M),nBlock,IDISK)
      END IF
#else
*     nTri=(nSize*(nSize+1))/2
      CALL DDAFILE(LUSBT,2,WORK(lg_M),nBlock,IDISK)
#endif


      END

      REAL*8 FUNCTION PSBMAT_FPRINT(lg_M,NM)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"

#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
#ifdef _GA_
        PSBMAT_FPRINT=SQRT(GA_DDOT(lg_M,lg_M))
#else
        MYRANK=GA_NODEID()
        NPROCS=GA_NNODES()
        ! get local stripes of RHS vectors
        CALL GA_Distribution (lg_M,myRank,iLo,iHi,jLo,jHi)
        DOTP=0.0D0
        IF (iLo.NE.0) THEN
          CALL GA_Access (lg_M,iLo,iHi,jLo,jHi,mM,LDM)
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          DOTP=DDOT_(NROW*NCOL,DBL_MB(mM),1,DBL_MB(mM),1)
          CALL GA_Release (lg_M,iLo,iHi,jLo,jHi)
        END IF
        CALL GADSUM_SCAL(DOTP)
        PSBMAT_FPRINT=SQRT(DOTP)
#endif
      ELSE
        nTri=(NM*(NM+1))/2
        PSBMAT_FPRINT=DNRM2_(nTri,WORK(lg_M),1)
      END IF
#else
      nTri=(NM*(NM+1))/2
      PSBMAT_FPRINT=DNRM2_(nTri,WORK(lg_M),1)
#endif
      END
