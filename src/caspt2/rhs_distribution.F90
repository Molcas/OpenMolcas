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
