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

      SUBROUTINE RHS_GET (NAS,NIS,lg_W,W)
      use definitions, only: iwp, wp
!SVC: this routine copies a global array to a local buffer
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
! SVC: when the _total_ size of a message exceeds 2**31-1 _bytes_,
! some implementations (e.g. MPICH, and thus also Intel MPI) fail.
! If this is the case, chop up the largest dimension and perform the
! GA_Get in batches smaller than 2**31-1 bytes (I took 2**30).
        MAX_MESG_SIZE = 2**27
        IF (NAS*NIS.GT.MAX_MESG_SIZE) THEN
          NIS_BATCH = MAX_MESG_SIZE / NAS
          IF (NIS_BATCH.EQ.0) THEN
            WRITE(u6,'(1X,A)') 'RHS_GET: NAS exceeds MAX_MESG_SIZE:'
            WRITE(u6,'(1X,I12,A,I12)') NAS, ' > ', MAX_MESG_SIZE
            CALL AbEnd()
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
