!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

      SUBROUTINE RHS_STRANS(NAS,NIS,ALPHA,lg_V1,lg_V2,ICASE,ISYM)
!SVC: this routine transforms RHS array V1 by multiplying on the left
!     with the S matrix and adds the result in V2: V2 <- V2 + alpha S*V1
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use fake_GA, only: GA_Arrays
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use constants, only: One
      use definitions, only: u6
#endif
      use definitions, only: iwp, wp
      IMPLICIT None
      integer(kind=iwp), intent(in):: NAS,NIS,lg_V1,lg_V2,ICASE,ISYM
      real(kind=wp), intent(in):: ALPHA

      real(kind=wp), Allocatable:: S(:)
      integer(kind=iwp) IDS, NS

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL(kind=iwp) bStat
      integer(kind=iwp) lg_S,myRank,iLoV1,iHiV1,jLoV1,jHiV1,            &
     &                              iLoV2,iHiV2,jLoV2,jHiV2,            &
     &                  NROW1,NROW2,NCOL1,NCOL2,mV1,LDV1,mV2,LDV2

      IF (Is_Real_Par()) THEN
        IF (ICASE.EQ.1 .OR. ICASE.EQ.4) THEN
!-SVC: if case is A or C, the S/ST matrices are loaded as global arrays,
!      then use the dgemm from GA to operate.
          CALL PSBMAT_GETMEM('S',lg_S,NAS)
          CALL PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
          CALL GA_DGEMM ('N','N',NAS,NIS,NAS,                           &
     &                   ALPHA,lg_S,lg_V1,One,lg_V2)
          bStat = GA_Destroy(lg_S)
        ELSE
!-SVC: if case is not A or C, the S/ST matrices are stored in replicate
!      fashion, and the RHS are stored as vertical stripes, so use
!      trimul on local memory, after accessing the local patch of the
!      vector.
          NS=(NAS*(NAS+1))/2
          CALL mma_allocate(S,NS,Label='S')
          IDS=IDSMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,S,NS,IDS)
!-SVC: get the local vertical stripes of the V1 and V2 vectors
          CALL GA_Sync()
          myRank = GA_NodeID()
          CALL GA_Distribution (lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
          CALL GA_Distribution (lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
          IF (jLoV1.NE.0.AND.jLoV2.NE.0) THEN
            NROW1=iHiV1-iLoV1+1
            NROW2=iHiV2-iLoV2+1
            NCOL1=jHiV1-jLoV1+1
            NCOL2=jHiV2-jLoV2+1
            IF (NCOL1.NE.NCOL2 .OR. NROW1.NE.NROW2 .OR.                 &
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
            CALL TRIMUL(NAS,NCOL1,ALPHA,S,                              &
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
        CALL TRIMUL(NAS,NIS,ALPHA,S,                                    &
     &              GA_Arrays(lg_V1)%A,NAS,                             &
     &              GA_Arrays(lg_V2)%A,NAS)
        CALL mma_deallocate(S)
#ifdef _MOLCAS_MPP_
      END IF
#include "macros.fh"
      unused_var(bStat)
#endif

      END SUBROUTINE RHS_STRANS
