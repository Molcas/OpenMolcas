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

      SUBROUTINE RHS_ACCESS (NAS,NIS,lg_W,iLo,iHi,jLo,jHi,MW)
!SVC: this routine gives a pointer to the process-local part of the RHS
!     If there is no valid local block, then the routine returns 0 for
!     iLo and jLo, and -1 for iHi and jHi. This way, loops from lower
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
! get the superindex ranges of this process's block
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)

        IF (iLo.NE.0 .AND. (iHi-iLo+1).NE.NAS) THEN
          WRITE(6,*) 'RHS_ACCESS: mismatch in range of the superindices'
          CALL AbEnd()
        END IF

! if the block is non-empty, get access to the block
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
