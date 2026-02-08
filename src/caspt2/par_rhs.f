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
* Copyright (C) 2011, Steven Vancoillie                                *
************************************************************************
************************************************************************
* Written by Steven Vancoillie, May 2011
* A set of subroutines that can handle RHS arrays in either a serial or
* parallel environment, depending on the situation.
************************************************************************
* --> when running serially, the RHS arrays are stored on LUSOLV and are
* loaded into the WORK array when needed.
* --> when running in parallel, the RHS arrays are stored on disk as
* disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
* number referring to the case, symmetry, and RHS vector respectively,
* and are loaded onto a global array when needed.
************************************************************************

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_INIT()
      use definitions, only: iwp, wp
      use caspt2_global, only: LURHS
      use caspt2_module, only: nSym, nISup, nASup, iOffRHS
      IMPLICIT None
      REAL(kind=wp) DUMMY(1)
      integer(kind=iwp) iDisk, iCase, iSym, NAS, NIS, NW, NRHS, iLo,
     &                  iHi, jLo, jHi

C-SVC: loop over symmetry/cases, get local patch of RHS, write, and then
C update the disk address in IOFFRHS
      IDISK=0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          IOFFRHS(ISYM,ICASE)=IDISK

          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NW=NAS*NIS

          IF (NW.EQ.0) CYCLE

          CALL RHS_DISTRIBUTION(NAS,NIS,iLo,iHi,jLo,jHi)
          NRHS=NAS*(jHi-jLo+1)
          CALL DDAFILE(LURHS(1),0,DUMMY,NRHS,IDISK)

        END DO
      END DO

      END SUBROUTINE RHS_INIT

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_FPRINT(CTYPE,IVEC)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero
      use caspt2_module, only: NSYM, NASUP, NINDEP, NISUP
      IMPLICIT None

      integer(kind=iwp), Intent(in):: IVEC
      CHARACTER(LEN=*), intent(in):: CTYPE

      real(kind=wp) :: FP(8)
      integer(kind=iwp) NROW, ICASE, ISYM, NAS, NIN, NIS, lg_W
      real(kind=wp), external:: RHS_DDOT

C-SVC: print out DNRM2 of the all RHS components
      NROW=0 ! dummy initialize
      DO ICASE=1,13
        DO ISYM=1,NSYM

          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)

          IF (CTYPE.EQ.'C') THEN
            NROW=NAS
          ELSE IF (CTYPE.EQ.'SR') THEN
            NROW=NIN
          ELSE
            WRITE(u6,'(1X,A)') 'RHS_FPRINT: invalid type: '//CTYPE
            CALL ABEND()
          END IF

          IF (NAS.NE.0 .AND. NIN.NE.0 .AND. NIS.NE.0) THEN
            CALL RHS_ALLO(NROW,NIS,lg_W)
            CALL RHS_READ(NROW,NIS,lg_W,iCASE,iSYM,iVEC)
            FP(ISYM)=SQRT(RHS_DDOT(NROW,NIS,lg_W,lg_W))
            CALL RHS_FREE(lg_W)
          ELSE
            FP(ISYM)=Zero
          END IF
        END DO
        WRITE(u6,'(1X,I2,1X,8F21.14)') ICASE, FP(1:NSYM)
      END DO

      END SUBROUTINE RHS_FPRINT

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_ZERO(IVEC)
      use definitions, only: iwp
      use constants, only: Zero
      use caspt2_module, only: NSYM,NASUP,NISUP
      IMPLICIT None
      Integer(kind=iwp), Intent(in):: IVEC

      Integer(kind=iwp) ICASE, ISYM, NAS, NIS, NW, lg_W

C-SVC: zero out the entire RHS vector on IVEC
      DO ICASE=1,13
        DO ISYM=1,NSYM

          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NW=NAS*NIS

          IF (NW.NE.0) THEN
            CALL RHS_ALLO(NAS,NIS,lg_W)
            CALL RHS_SCAL(NAS,NIS,lg_W,Zero)
            CALL RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
            CALL RHS_FREE(lg_W)
          END IF
        END DO
      END DO

      END SUBROUTINE RHS_ZERO

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_ALLO (NAS,NIS,lg_W)
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: Allocate_GA_Array
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS,NIS
      integer(kind=iwp), intent(out):: lg_W

      integer(kind=iwp) NW
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"

      IF (Is_Real_Par()) THEN
        CALL GA_CREATE_STRIPED ('V',NAS,NIS,'RHS',LG_W)
      ELSE
#endif
        NW=NAS*NIS
        lg_W=Allocate_GA_Array(NW,'RHS')
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_ALLO

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_FREE (lg_W)
CSVC: this routine writes the RHS array to disk
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: Deallocate_GA_Array
      IMPLICIT None
      integer(kind=iwp), intent(inout):: lg_w
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL(kind=iwp) bStat

      IF (Is_Real_Par()) THEN
CSVC: Destroy the global array
        bStat=GA_Destroy(lg_W)
      ELSE
#endif
        Call Deallocate_GA_Array(lg_W)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_FREE

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_DISTRIBUTION (NAS,NIS,iLo,iHi,jLo,jHi)
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT None
      integer(kind=iwp), intent(in) :: NAS,NIS
      integer(kind=iwp), intent(out) :: iLo,iHi,jLo,jHi


#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) MYRANK, NPROCS, NBASE, NREST
#endif

      iLo=1
      iHi=NAS

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        MYRANK=GA_NODEID()
        NPROCS=GA_NNODES()
        NBASE=NIS/NPROCS
        NREST=NIS-NBASE*NPROCS
        IF (MYRANK.LT.NREST) THEN
          jLo=MYRANK*(NBASE+1)+1
          jHi=jLo+NBASE
        ELSE
          jLo=NREST*(NBASE+1)+(MYRANK-NREST)*NBASE+1
          jHi=jLo+NBASE-1
        END IF
      ELSE
#endif
        jLo=1
        jHi=NIS
#ifdef _MOLCAS_MPP_
      END IF
#endif
      END SUBROUTINE RHS_DISTRIBUTION

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_ACCESS (NAS,NIS,lg_W,iLo,iHi,jLo,jHi,MW)
CSVC: this routine gives a pointer to the process-local part of the RHS
C     If there is no valid local block, then the routine returns 0 for
C     iLo and jLo, and -1 for iHi and jHi. This way, loops from lower
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS,NIS,lg_W
      integer(kind=iwp), intent(out):: iLo,iHi,jLo,jHi,MW
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) myRank, LDW
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
* get the superindex ranges of this process's block
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)

        IF (iLo.NE.0 .AND. (iHi-iLo+1).NE.NAS) THEN
          WRITE(6,*) 'RHS_ACCESS: mismatch in range of the superindices'
          CALL AbEnd()
        END IF

* if the block is non-empty, get access to the block
        IF (iLo.GT.0 .AND. jLo.GT.0) THEN
          IF (iHi-iLo+1.NE.NAS) THEN
            WRITE(6,*) 'RHS_ACCESS: Error: NAS mismatch, abort...'
            CALL ABEND()
          END IF
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,MW,LDW)
          IF (LDW.NE.NAS) THEN
            WRITE(6,*) 'RHS_ACCESS: assert NAS=LDW failed, abort'
            CALL AbEnd()
          END IF
        END IF
      ELSE
#endif
        iLo=1
        iHi=NAS
        jLo=1
        jHi=NIS
        MW=lg_W
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_ACCESS
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_RELEASE (lg_W,iLo,iHi,jLo,jHi)
CSVC: this routine releases a local block back to the global array
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT None
      integer(kind=iwp), Intent(inout):: lg_W,iLo,iHi,jLo,jHi
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"

      IF (Is_Real_Par()) THEN
        IF (iLo.GT.0 .AND. jLo.GT.0) THEN
          CALL GA_Release (lg_W,iLo,iHi,jLo,jHi)
        END IF
      END IF
#else
#include "macros.fh"
C Avoid unused argument warnings
      unused_var(lg_W)
      unused_var(iLo)
      unused_var(iHi)
      unused_var(jLo)
      unused_var(jHi)
#endif

      END SUBROUTINE RHS_RELEASE
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_RELEASE_UPDATE (lg_W,iLo,iHi,jLo,jHi)
CSVC: this routine releases a local block that was written to back to
C the global array
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT None
      integer(kind=iwp), Intent(inout):: lg_W,iLo,iHi,jLo,jHi
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"

      IF (Is_Real_Par()) THEN
        IF (iLo.GT.0 .AND. jLo.GT.0) THEN
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
      END IF
#else
#include "macros.fh"
C Avoid unused argument warnings
      unused_var(lg_W)
      unused_var(iLo)
      unused_var(iHi)
      unused_var(jLo)
      unused_var(jHi)
#endif

      END SUBROUTINE RHS_RELEASE_UPDATE
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_GET (NAS,NIS,lg_W,W)
      use definitions, only: iwp, wp
CSVC: this routine copies a global array to a local buffer
#ifdef _MOLCAS_MPP_
      use definitions, only: u6
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      integer(kind=iwp), Intent(In):: NAS,NIS,lg_W
      real(kind=wp), Intent(Out):: W(NAS*NIS)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) MAX_MESG_SIZE, NIS_BATCH, NIS_STA, NIS_END, IOFF

      IF (Is_Real_Par()) THEN
C SVC: when the _total_ size of a message exceeds 2**31-1 _bytes_,
C some implementations (e.g. MPICH, and thus also Intel MPI) fail.
C If this is the case, chop up the largest dimension and perform the
C GA_Get in batches smaller than 2**31-1 bytes (I took 2**30).
        MAX_MESG_SIZE = 2**27
        IF (NAS*NIS.GT.MAX_MESG_SIZE) THEN
          NIS_BATCH = MAX_MESG_SIZE / NAS
          IF (NIS_BATCH.EQ.0) THEN
            WRITE(u6,'(1X,A)') 'RHS_GET: NAS exceeds MAX_MESG_SIZE:'
            WRITE(u6,'(1X,I12,A,I12)') NAS, ' > ', MAX_MESG_SIZE
            CALL AbEnd
          END IF
          DO NIS_STA=1,NIS,NIS_BATCH
            NIS_END=MIN(NIS_STA+NIS_BATCH-1,NIS)
            IOFF=NAS*(NIS_STA-1)+1
            CALL GA_Get (lg_W,1,NAS,NIS_STA,NIS_END,W(IOFF),NAS)
          END DO
        ELSE
          CALL GA_Get (lg_W,1,NAS,1,NIS,W,NAS)
        END IF
      ELSE
#endif
        CALL DCOPY_(NAS*NIS,GA_Arrays(lg_W)%A,1,W,1)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_GET
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_PUT (NAS,NIS,lg_W,W)
CSVC: this routine copies a local buffer to a global array
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
      use definitions, only: u6
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      integer(kind=iwp), Intent(In):: NAS,NIS,lg_W
      real(kind=wp), Intent(in):: W(NAS*NIS)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) MAX_MESG_SIZE, NIS_BATCH, NIS_STA, NIS_END, IOFF

      IF (Is_Real_Par()) THEN
        IF (KING()) THEN
C SVC: when the _total_ size of a message exceeds 2**31-1 _bytes_,
C some implementations (e.g. MPICH, and thus also Intel MPI) fail.
C If this is the case, chop up the largest dimension and perform the
C GA_Put in batches smaller than 2**31-1 bytes (I took 2**27 elements
C which is 2**30 bytes).
          MAX_MESG_SIZE = 2**27
          IF (NAS*NIS.GT.MAX_MESG_SIZE) THEN
            NIS_BATCH = MAX_MESG_SIZE / NAS
            IF (NIS_BATCH.EQ.0) THEN
              WRITE(u6,'(1X,A)') 'RHS_GET: NAS exceeds MAX_MESG_SIZE:'
              WRITE(u6,'(1X,I12,A,I12)') NAS, ' > ', MAX_MESG_SIZE
              CALL AbEnd
            END IF
            DO NIS_STA=1,NIS,NIS_BATCH
              NIS_END=MIN(NIS_STA+NIS_BATCH-1,NIS)
              IOFF=NAS*(NIS_STA-1)+1
              CALL GA_Put (lg_W,1,NAS,NIS_STA,NIS_END,W(IOFF),NAS)
            END DO
          ELSE
            CALL GA_Put (lg_W,1,NAS,1,NIS,W,NAS)
          END IF
        END IF
      ELSE
#endif
        CALL DCOPY_(NAS*NIS,W,1,GA_Arrays(lg_W)%A,1)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_PUT

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_ADD (NAS,NIS,lg_W,W)
CSVC: this routine adds to the local part of a global RHS array the
      use definitions, only: iwp, wp
      use constants, only: One
Cmatching part of a replicate array.
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      integer(kind=iwp), Intent(in):: NAS,NIS,lg_W
      real(kind=wp), Intent(In):: W(NAS,*)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) myRank,iLo,iHi,jLo,jHi,NW,mW,LDW

      IF (Is_Real_Par()) THEN
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NW=(iHi-iLo+1)*(jHi-jLo+1)
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,mW,LDW)
          CALL DAXPY_(NW,One,W(iLo,jLo),1,DBL_MB(mW),1)
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
      ELSE
#endif
        CALL DAXPY_(NAS*NIS,One,W,1,GA_Arrays(lg_W)%A,1)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_ADD

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_READ_C (lg_W,iCASE,iSYM,iVEC)
      use definitions, only: iwp
      use caspt2_module, only: NASUP, NISUP
      IMPLICIT None
      integer(kind=iwp), Intent(In):: lg_W,iCASE,iSYM,iVEC
      integer(kind=iwp) NAS,NIS
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      CALL RHS_READ (NAS,NIS,lg_W,ICASE,ISYM,IVEC)
      END SUBROUTINE RHS_READ_C

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_READ_SR (lg_W,iCASE,iSYM,iVEC)
      use definitions, only: iwp
      use caspt2_module, only: NINDEP, NISUP
      IMPLICIT None
      integer(kind=iwp), Intent(In):: lg_W,iCASE,iSYM,iVEC
      integer(kind=iwp) NIN,NIS
      NIN=NINDEP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      CALL RHS_READ (NIN,NIS,lg_W,ICASE,ISYM,IVEC)
      END SUBROUTINE RHS_READ_SR

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_READ (NIN,NIS,lg_W,iCASE,iSYM,iVEC)
CSVC: this routine reads an RHS array in SR format from disk
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_global, only: LURHS
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: IOFFRHS
      IMPLICIT None
      integer(kind=iwp), Intent(In):: NIN,NIS,lg_W,iCASE,iSYM,iVEC

      integer(kind=iwp) IDISK, NWPROC
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) myRank,ISTA,IEND,JSTA,JEND,mpt_W,LDW

      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_W,myRank,ISTA,IEND,JSTA,JEND)
        IF (IEND-ISTA+1.EQ.NIN .AND. ISTA.GT.0) THEN
          CALL GA_Access (lg_W,ISTA,IEND,JSTA,JEND,mpt_W,LDW)
          IF (LDW.NE.NIN) THEN
            WRITE(6,*) 'RHS_READ: Assumption NIN==LDW wrong'
            CALL AbEnd()
          END IF
          NWPROC=NIN*(JEND-JSTA+1)
          IDISK=IOFFRHS(ISYM,ICASE)
          CALL DDAFILE(LURHS(IVEC),2,DBL_MB(mpt_W),NWPROC,IDISK)
          CALL GA_Release_Update (lg_W,ISTA,IEND,JSTA,JEND)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        NWPROC=NIN*NIS
        IDISK=IOFFRHS(ISYM,ICASE)
        CALL DDAFILE(LURHS(IVEC),2,GA_Arrays(lg_W)%A,NWPROC,IDISK)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_READ

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_SAVE_C (lg_W,iCASE,iSYM,iVEC)
      use definitions, only: iwp
      use caspt2_module, only: NASUP, NISUP
      IMPLICIT None
      integer(kind=iwp), intent(in):: lg_W,iCASE,iSYM,iVEC
      integer(kind=iwp) NAS, NIS
      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      CALL RHS_SAVE (NAS,NIS,lg_W,ICASE,ISYM,IVEC)
      END SUBROUTINE RHS_SAVE_C

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_SAVE_SR (lg_W,iCASE,iSYM,iVEC)
      use definitions, only: iwp
      use caspt2_module, only: NINDEP, NISUP
      IMPLICIT None
      integer(kind=iwp), intent(in):: lg_W,iCASE,iSYM,iVEC
      integer(kind=iwp) NIN, NIS
      NIN=NINDEP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)
      CALL RHS_SAVE (NIN,NIS,lg_W,ICASE,ISYM,IVEC)
      END SUBROUTINE RHS_SAVE_SR

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_SAVE (NIN,NIS,lg_W,iCASE,iSYM,iVEC)
      use definitions, only: iwp
CSVC: this routine reads an RHS array in SR format from disk
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use definitions, only: u6
#endif
      use caspt2_global, only: LURHS
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: IOFFRHS
      IMPLICIT None
      integer(kind=iwp), intent(in):: NIN,NIS,lg_W,iCASE,iSYM,iVEC

      integer(kind=iwp) IDISK, NW
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) myRank,ISTA,IEND,JSTA,JEND,mpt_W,LDW,NWPROC
#include "global.fh"
#include "mafdecls.fh"

      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_W,myRank,ISTA,IEND,JSTA,JEND)
        IF (IEND-ISTA+1.EQ.NIN .AND. ISTA.GT.0) THEN
          CALL GA_Access (lg_W,ISTA,IEND,JSTA,JEND,mpt_W,LDW)
          IF (LDW.NE.NIN) THEN
            WRITE(u6,*) 'RHS_SAVE: Assumption NIN==LDW wrong'
            CALL AbEnd()
          END IF
          NWPROC=NIN*(JEND-JSTA+1)
          IDISK=IOFFRHS(ISYM,ICASE)
          CALL DDAFILE(LURHS(IVEC),1,DBL_MB(mpt_W),NWPROC,IDISK)
          CALL GA_Release (lg_W,ISTA,IEND,JSTA,JEND)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        NW=NIN*NIS
        IDISK=IOFFRHS(ISYM,ICASE)
        CALL DDAFILE(LURHS(IVEC),1,GA_Arrays(lg_W)%A,NW,IDISK)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SAVE

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_SCATTER (LDW,lg_W,Buff,idxW,nBuff)
CSVC: this routine scatters + adds values of a buffer array into the RHS
C     array at positions given by the buffer index array.
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use stdalloc, only: mma_allocate, mma_deallocate
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      Integer(kind=iwp), intent(in):: LDW,lg_W,nBuff
      real(kind=wp), intent(in):: Buff(nBuff)
      Integer(kind=iwp), intent(in):: idxW(nBuff)

      Integer(kind=iwp) I
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      Integer(kind=iwp), allocatable:: TMPW1(:), TMPW2(:)
#else
#include "macros.fh"
      unused_var(LDW)
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
CSVC: global array RHS matrix expects 2 index buffers
        CALL mma_allocate(TMPW1,nBuff,Label='TMPW1')
        CALL mma_allocate(TMPW2,nBuff,Label='TMPW2')
        DO I=1,nBuff
          TMPW2(I)=(idxW(I)-1)/LDW+1
          TMPW1(I)=idxW(I)-LDW*(TMPW2(I)-1)
        END DO
        CALL GA_Scatter_Acc (lg_W,Buff,TMPW1,TMPW2,nBuff,1.0D0)
        CALL mma_deallocate(TMPW1)
        CALL mma_deallocate(TMPW2)
      ELSE
#endif
        DO I=1,nBuff
          GA_Arrays(lg_W)%A(idxW(I)) =
     &      GA_Arrays(lg_W)%A(idxW(I)) + BUFF(I)
        END DO
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SCATTER

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE DRA2SOLV (NAS,NIS,iCASE,iSYM,iVEC)
CSVC: FIXME: this temporary routine copies the RHS arrays from DRAs to
C     LUSOLV and should be removed once the full parallelization is in
C     place and transition is no longer needed.
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
      use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
      use definitions, only: wp, u6
#endif
      use caspt2_global, only: IDSCT, LUSOLV
      use EQSOLV, only: MXSCT
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: MXCASE
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS,NIS,iCASE,iSYM,iVEC

      integer(kind=iwp) IDISK, lg_W
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
*     LOGICAL bStat
      real(kind=wp), allocatable:: TMPW(:)
      integer(kind=iwp) iMax, NCOL, NW, ISTA, IEND
#endif

CSVC: Read the global array from disk
      CALL RHS_ALLO (NAS,NIS,lg_W)
      CALL RHS_READ (NAS,NIS,lg_W,ICASE,ISYM,IVEC)

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
CSVC: Only the master process writes to LUSOLV!!
CSVC: be careful to only call one-sided operations
        IF (KING()) THEN
CSVC: write the global array to LUSOLV
C     when all routines have been adapted, this can go into SYNRHS and
C     later it should be completely removed when VCUTIL and SGM are
C     adapted for handling global + disk resident arrays then, also
C     remove iCASE,iSYM,iVEC from the call, as they are no longer needed
C     in that case.
          CALL mma_MaxDBLE(iMax)
C-SVC: GA_Get does not like large buffer sizes, put upper limit at 1GB
          iMax=MIN(NINT(0.95D0*iMax),134217728)
          NCOL=MIN(iMAX,NAS*NIS)/NAS
          IF (NCOL.LE.0) THEN
            WRITE(u6,*) 'Not enough memory in DRA2SOLV, aborting...'
            CALL AbEnd()
          END IF
          NW=NAS*NCOL
          CALL mma_allocate(TMPW,NW,LABEL='TMPW')
CSVC: Write local array to LUSOLV
          IDISK=IDSCT(1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
          DO ISTA=1,NIS,NCOL
            IEND=MIN(ISTA+NCOL-1,NIS)
            CALL GA_Get (lg_W,1,NAS,ISTA,IEND,TMPW,NAS)
            CALL DDAFILE(LUSOLV,1,TMPW,NAS*(IEND-ISTA+1),IDISK)
          END DO
          CALL mma_deallocate(TMPW)
        END IF
        CALL GASync()
CSVC: Destroy the global array
*       bStat=GA_Destroy(lg_W)
      ELSE
#endif
      IDISK=IDSCT(1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
      CALL DDAFILE(LUSOLV,1,GA_Arrays(lg_W)%A,NAS*NIS,IDISK)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      CALL RHS_FREE (lg_W)

      END SUBROUTINE DRA2SOLV

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE SOLV2DRA (NAS,NIS,iCASE,iSYM,iVEC)
CSVC: FIXME: this temporary routine copies the RHS arrays from DRAs to
C     LUSOLV and should be removed once the full parallelization is in
C     place and transition is no longer needed.
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      use definitions, only: wp, u6
      USE Para_Info, ONLY: Is_Real_Par, King
      use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
#endif
      use caspt2_global, only: LUSOLV, IDSCT
      use EQSOLV, only: MXSCT
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: MXCASE
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS,NIS,iCASE,iSYM,iVEC

      integer(kind=iwp) IDISK, lg_W, NW
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
*     LOGICAL bStat
      real(kind=wp), allocatable:: TMPW(:)
      integer(kind=iwp) iMax, NCOL, ISTA, IEND
#endif

      CALL RHS_ALLO (NAS,NIS,lg_W)

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
CSVC: Read the global array from disk
CSVC: Only the master process writes to LUSOLV!!
CSVC: be careful to only call one-sided operations
        IF (KING()) THEN
CSVC: write the LUSOLV array to global RHS array
C     later it should be completely removed when everything is parallel
          CALL mma_MaxDBLE(iMax)
C-SVC: GA_Get does not like large buffer sizes, put upper limit at 1GB
          iMax=MIN(NINT(0.95D0*iMax),134217728)
          NCOL=MIN(iMAX,NAS*NIS)/NAS
          IF (NCOL.LE.0) THEN
            WRITE(u6,*) 'Not enough memory in SOLV2DRA, aborting...'
            CALL AbEnd()
          END IF
          NW=NAS*NCOL
          CALL mma_allocate(TMPW,NW,Label='TMPW')
CSVC: Read local array from LUSOLV
          IDISK=IDSCT(1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
          DO ISTA=1,NIS,NCOL
            IEND=MIN(ISTA+NCOL-1,NIS)
            CALL DDAFILE(LUSOLV,2,TMPW,NAS*(IEND-ISTA+1),IDISK)
            CALL GA_Put (lg_W,1,NAS,ISTA,IEND,TMPW,NAS)
          END DO
          CALL mma_deallocate(TMPW)
        END IF
        CALL GASync()
CSVC: Destroy the global array
*       bStat=GA_Destroy(lg_W)
      ELSE
#endif
        NW=NAS*NIS
        IDISK=IDSCT(1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
        CALL DDAFILE(LUSOLV,2,GA_Arrays(lg_W)%A,NW,IDISK)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      CALL RHS_SAVE (NAS,NIS,lg_W,ICASE,ISYM,IVEC)
      CALL RHS_FREE (lg_W)

      END SUBROUTINE SOLV2DRA

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_SCAL (NAS,NIS,lg_W,FACT)
CSVC: this routine multiplies the RHS array with FACT
      use definitions, only: iwp, wp
      use constants, only: Zero, One
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS,NIS,lg_W
      real(kind=wp), intent(in):: FACT
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        IF (FACT.EQ.Zero) THEN
C          CALL GA_Fill (lg_W,0.0D0)
           CALL GA_Zero (lg_W)
        ELSE
          IF (FACT.NE.One) THEN
            CALL GA_Scale (lg_W,FACT)
          END IF
        END IF
      ELSE
#endif
        IF(FACT.EQ.Zero) THEN
            CALL DCOPY_(NAS*NIS,[Zero],0,GA_Arrays(lg_W)%A,1)
        ELSE
          IF(FACT.NE.One) THEN
            CALL DSCAL_(NAS*NIS,FACT,GA_Arrays(lg_W)%A,1)
          END IF
        END IF
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SCAL

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_SR2C (ITYP,IREV,NAS,NIS,NIN,lg_V1,lg_V2,
     &                     ICASE,ISYM)
      use definitions, only: iwp, wp, u6
      use Constants, only: Zero, One
CSVC: this routine transforms the RHS arrays from SR format (V1) to C
C     format (V2) (IREV=0) and back (IREV=1), with ITYP specifying if
C     only the T matrix is used (ITYP=0) or the product of S and T
C     (ITYP=1).
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDTMAT, IDSTMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      integer(kind=iwp), intent(in)::ITYP,IREV,NAS,NIS,NIN,lg_V1,lg_V2,
     &                               ICASE,ISYM

      Real(kind=wp), Allocatable:: T(:)
      integer(kind=iwp) IDT

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL(kind=iwp) bStat
      integer(kind=iwp) lg_T, myRank,iLoV1,iHiV1,jLoV1,jHiV1,
     &                               iLoV2,iHiV2,jLoV2,jHiV2,
     &                  NROW1, NROW2, NCOL1, NCOL2,
     &                  mV1,LDV1,mV2,LDV2

      IF (Is_Real_Par()) THEN
        IF (ICASE.EQ.1 .OR. ICASE.EQ.4) THEN
C-SVC: if case is A or C, the S/ST matrices are loaded as global arrays,
C      then use the dgemm from GA to operate.
          CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
          IF (ITYP.EQ.0) THEN
            CALL PSBMAT_READ ('T',iCase,iSym,lg_T,NAS*NIN)
          ELSE IF (ITYP.EQ.1) THEN
            CALL PSBMAT_READ ('M',iCase,iSym,lg_T,NAS*NIN)
          ELSE
            WRITE(u6,*) 'RHS_SR2C: invalid type = ', ITYP
            CALL AbEnd()
          END IF

          IF (IREV.EQ.0) THEN
            CALL GA_DGEMM ('N','N',NAS,NIS,NIN,
     &                     One,lg_T,lg_V1,Zero,lg_V2)
          ELSE
            CALL GA_DGEMM ('T','N',NIN,NIS,NAS,
     &                     One,lg_T,lg_V2,Zero,lg_V1)
          END IF
          bStat = GA_Destroy(lg_T)
        ELSE
C-SVC: if case is not A or C, the S/ST matrices are stored in replicate
C      fashion, and the RHS are stored as vertical stripes, so use dgemm
C      on local memory, after accessing the local patch of the vector.
          CALL mma_allocate(T,NAS*NIN,Label='T')
          IF (ITYP.EQ.0) THEN
            IDT=IDTMAT(ISYM,ICASE)
          ELSE IF (ITYP.EQ.1) THEN
            IDT=IDSTMAT(ISYM,ICASE)
          ELSE
            WRITE(u6,*) 'RHS_SR2C: invalid type = ', ITYP
            CALL AbEnd()
          END IF
          CALL DDAFILE(LUSBT,2,T,NAS*NIN,IDT)
C-SVC: get the local vertical stripes of the V1 and V2 vectors
          CALL GA_Sync()
          myRank = GA_NodeID()
          CALL GA_Distribution (lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
          CALL GA_Distribution (lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
          IF (jLoV1.NE.0.AND.jLoV2.NE.0) THEN
            NROW1=iHiV1-iLoV1+1
            NROW2=iHiV2-iLoV2+1
            NCOL1=jHiV1-jLoV1+1
            NCOL2=jHiV2-jLoV2+1
            IF (NCOL1.NE.NCOL2 .OR. NROW1.NE.NIN .OR. NROW2.NE.NAS) THEN
              WRITE(u6,*) 'RHS_SR2C: inconsistent stripe size'
              WRITE(u6,'(A,I3)') 'ICASE = ', ICASE
              WRITE(u6,'(A,I3)') 'ISYM  = ', ISYM
              WRITE(u6,'(A,2I6)') 'NCOL1, NCOL2 = ', NCOL1, NCOL2
              WRITE(u6,'(A,2I6)') 'NROW1, NIN   = ', NROW1, NIN
              WRITE(u6,'(A,2I6)') 'NROW2, NAS   = ', NROW2, NAS
              CALL AbEnd()
            END IF
            CALL GA_Access (lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
            CALL GA_Access (lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
            IF (IREV.EQ.0) THEN
              CALL DGEMM_('N','N',NAS,NCOL1,NIN,
     &                    One,T,NAS,DBL_MB(mV1),LDV1,
     &                    Zero,DBL_MB(mV2),LDV2)
            ELSE
              CALL DGEMM_('T','N',NIN,NCOL1,NAS,
     &                    One,T,NAS,DBL_MB(mV2),LDV2,
     &                    Zero,DBL_MB(mV1),LDV1)
*             WRITE(6,*) 'Fingerprint =', RHS_DDOT(NAS,NIN,lg_V1,lg_V1)
            END IF
            CALL GA_Release_Update (lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
            CALL GA_Release_Update (lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
          END IF
          CALL mma_deallocate(T)
          CALL GA_Sync()
        END IF
      ELSE
#endif
        CALL mma_allocate(T,NAS*NIN,Label='T')
        IF (ITYP.EQ.0) THEN
          IDT=IDTMAT(ISYM,ICASE)
        ELSE IF (ITYP.EQ.1) THEN
          IDT=IDSTMAT(ISYM,ICASE)
        ELSE
          WRITE(u6,*) 'RHS_SR2C: invalid type = ', ITYP
          CALL AbEnd()
        END IF
        CALL DDAFILE(LUSBT,2,T,NAS*NIN,IDT)
        IF (IREV.EQ.0) THEN
          CALL DGEMM_('N','N',NAS,NIS,NIN,
     &                One,T,NAS,GA_Arrays(lg_V1)%A,NIN,
     &                Zero,GA_Arrays(lg_V2)%A,NAS)
        ELSE
          CALL DGEMM_('T','N',NIN,NIS,NAS,
     &                One,T,NAS,GA_Arrays(lg_V2)%A,NAS,
     &                Zero,GA_Arrays(lg_V1)%A,NIN)
        END IF
        CALL mma_deallocate(T)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SR2C

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_STRANS(NAS,NIS,ALPHA,lg_V1,lg_V2,ICASE,ISYM)
CSVC: this routine transforms RHS array V1 by multiplying on the left
C     with the S matrix and adds the result in V2: V2 <- V2 + alpha S*V1
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      use definitions, only: u6
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS,NIS,lg_V1,lg_V2,ICASE,ISYM
      real(kind=wp), intent(in):: ALPHA

      real(kind=wp), Allocatable:: S(:)
      integer(kind=iwp) IDS, NS

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL(kind=iwp) bStat
      integer(kind=iwp) lg_S,myRank,iLoV1,iHiV1,jLoV1,jHiV1,
     &                              iLoV2,iHiV2,jLoV2,jHiV2,
     &                  NROW1,NROW2,NCOL1,NCOL2,mV1,LDV1,mV2,LDV2

      IF (Is_Real_Par()) THEN
        IF (ICASE.EQ.1 .OR. ICASE.EQ.4) THEN
C-SVC: if case is A or C, the S/ST matrices are loaded as global arrays,
C      then use the dgemm from GA to operate.
          CALL PSBMAT_GETMEM('S',lg_S,NAS)
          CALL PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
          CALL GA_DGEMM ('N','N',NAS,NIS,NAS,
     &                   ALPHA,lg_S,lg_V1,1.0D0,lg_V2)
          bStat = GA_Destroy(lg_S)
        ELSE
C-SVC: if case is not A or C, the S/ST matrices are stored in replicate
C      fashion, and the RHS are stored as vertical stripes, so use
C      trimul on local memory, after accessing the local patch of the
C      vector.
          NS=(NAS*(NAS+1))/2
          CALL mma_allocate(S,NS,Label='S')
          IDS=IDSMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,S,NS,IDS)
C-SVC: get the local vertical stripes of the V1 and V2 vectors
          CALL GA_Sync()
          myRank = GA_NodeID()
          CALL GA_Distribution (lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
          CALL GA_Distribution (lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
          IF (jLoV1.NE.0.AND.jLoV2.NE.0) THEN
            NROW1=iHiV1-iLoV1+1
            NROW2=iHiV2-iLoV2+1
            NCOL1=jHiV1-jLoV1+1
            NCOL2=jHiV2-jLoV2+1
            IF (NCOL1.NE.NCOL2 .OR. NROW1.NE.NROW2 .OR.
     &          NROW1.NE.NAS) THEN
              WRITE(u6,*) 'RHS_STRANS: inconsistent stripe size'
              WRITE(u6,'(A,I3)') 'ICASE = ', ICASE
              WRITE(u6,'(A,I3)') 'ISYM  = ', ISYM
              WRITE(u6,'(A,2I6)') 'NCOL1, NCOL2 = ', NCOL1, NCOL2
              WRITE(u6,'(A,2I6)') 'NROW1, NROW2 = ', NROW1, NROW2
              CALL AbEnd()
            END IF
            CALL GA_Access (lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
            CALL GA_Access (lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
            CALL TRIMUL(NAS,NCOL1,ALPHA,S,
     &                  DBL_MB(mV1),LDV1,DBL_MB(mV2),LDV2)
            CALL GA_Release_Update (lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
            CALL GA_Release_Update (lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
          END IF
          CALL GA_Sync()
          CALL mma_deallocate(S)
        END IF
      ELSE
#endif
        NS=(NAS*(NAS+1))/2
        CALL mma_allocate(S,NS,Label='S')
        IDS=IDSMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,2,S,NS,IDS)
        CALL TRIMUL(NAS,NIS,ALPHA,S,
     &              GA_Arrays(lg_V1)%A,NAS,
     &              GA_Arrays(lg_V2)%A,NAS)
        CALL mma_deallocate(S)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_STRANS

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      FUNCTION RHS_DDOT(NAS,NIS,lg_V1,lg_V2)
CSVC: this routine computes the DDOT of the RHS arrays V1 and V2
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      real(kind=wp) RHS_DDOT
      Integer(kind=iwp), intent(in):: NAS,NIS,lg_V1,lg_V2
      real(kind=wp), external:: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"

      IF (Is_Real_Par()) THEN
        RHS_DDOT = GA_DDOT(lg_V1,lg_V2)
      ELSE
#endif
        RHS_DDOT = DDOT_(NAS*NIS,GA_Arrays(lg_V1)%A,1,
     &                           GA_Arrays(lg_V2)%A,1)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END FUNCTION RHS_DDOT

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_DAXPY (NAS,NIS,ALPHA,lg_V1,lg_V2)
CSVC: this routine computes product ALPHA * V1 and adds to V2
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS, NIS
      real(kind=wp), intent(in):: ALPHA
      integer(kind=iwp), intent(in):: lg_V1, lg_V2

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) myRank,iLoV1,iHiV1,jLoV1,jHiV1,
     &                         iLoV2,iHiV2,jLoV2,jHiV2,
     &                  NV1,NV2,mV1,LDV1,mV2,LDV2

      IF (Is_Real_Par()) THEN
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
        CALL GA_Distribution (lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
        IF (iLoV1.NE.0.AND.iLoV2.NE.0) THEN
          NV1=(iHiV1-iLoV1+1)*(jHiV1-jLoV1+1)
          NV2=(iHiV2-iLoV2+1)*(jHiV2-jLoV2+1)
          IF (NV1.NE.NV2) CALL AbEnd()
          CALL GA_Access (lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
          CALL GA_Access (lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
          ! V2 <- alpha*V1 + V2
          CALL DAXPY_(NV1,ALPHA,DBL_MB(mV1),1,DBL_MB(mV2),1)
          CALL GA_Release_Update (lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
          CALL GA_Release (lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
        END IF
      ELSE
#endif
        CALL DAXPY_(NAS*NIS,ALPHA,GA_Arrays(lg_V1)%A,1,
     &                            GA_Arrays(lg_V2)%A,1)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_DAXPY

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_RESDIA(NIN,NIS,lg_W,DIN,DIS,DOVL)
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use constants, only: Zero
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None

      integer(kind=iwp), intent(in):: NIN,NIS,lg_W
      real(kind=wp), Intent(in):: DIN(*),DIS(*)
      real(kind=wp), Intent(Out):: DOVL

C Apply the resolvent of the diagonal part of H0 to an RHS array

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) myRank,iLo,iHi,jLo,jHi,NROW,NCOL,mW,LDW

      IF (Is_Real_Par()) THEN
        DOVL=Zero
        CALL GA_Sync()
        myRank = GA_NodeID()
C-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,mW,LDW)
          CALL RESDIA(NROW,NCOL,DBL_MB(mW),LDW,DIN(iLo),DIS(jLo),DOVL)
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
        CALL GA_Sync()
        CALL GAdSUM_SCAL(DOVL)
      ELSE
#endif
        CALL RESDIA(NIN,NIS,GA_Arrays(lg_W)%A,NIN,DIN,DIS,DOVL)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_RESDIA

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_SGMDIA(NIN,NIS,lg_W,DIN,DIS)
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      IMPLICIT None

      integer(kind=iwp), intent(in):: NIN,NIS,lg_W
      real(kind=wp), Intent(in):: DIN(*),DIS(*)

C Apply the resolvent of the diagonal part of H0 to an RHS array

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer(kind=iwp) myRank,iLo,iHi,jLo,jHi,NROW,NCOL,mW,LDW

      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
C-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,mW,LDW)
          CALL SGMDIA(NROW,NCOL,DBL_MB(mW),LDW,DIN(iLo),DIS(jLo))
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        CALL SGMDIA(NIN,NIS,GA_Arrays(lg_W)%A,NIN,DIN,DIS)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SGMDIA
