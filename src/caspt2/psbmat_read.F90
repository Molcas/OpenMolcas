!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
! WRAPPERS FOR PARALLEL S AND B MATRIX ROUTINES
!
      SUBROUTINE PSBMAT_READ(cNAME,iCase,iSym,lg_M,nSize)
!SVC20100902: read the disk array stored as cName+iSym using DRA
! interface into global array lg_M, or if replicate or serial, read from
! LUSBT into WORK(lg_M)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use caspt2_global, only: LUH0T
#endif
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT, IDBMAT, IDTMAT, IDSTMAT
      use fake_ga, only: GA_arrays
      use definitions, only: iwp
      IMPLICIT None
      INTEGER(kind=iwp) iCASE,iSym,lg_M,nSize
      CHARACTER(LEN=*) cNAME

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      INTEGER(kind=iwp) LU, myRank, iSTA, IEND, JSTA, JEND, mpt_M, LDM
#endif
      INTEGER(kind=iwp) IDISK, nBlock


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
