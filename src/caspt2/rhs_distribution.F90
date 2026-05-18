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

subroutine RHS_DISTRIBUTION(NAS,NIS,iLo,iHi,jLo,jHi)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NAS, NIS
integer(kind=iwp), intent(out) :: iLo, iHi, jLo, jHi
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: MYRANK, NBASE, NPROCS, NREST
#include "global.fh"
#include "mafdecls.fh"
#endif

iLo = 1
iHi = NAS

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  MYRANK = GA_NODEID()
  NPROCS = GA_NNODES()
  NBASE = NIS/NPROCS
  NREST = NIS-NBASE*NPROCS
  if (MYRANK < NREST) then
    jLo = MYRANK*(NBASE+1)+1
    jHi = jLo+NBASE
  else
    jLo = NREST*(NBASE+1)+(MYRANK-NREST)*NBASE+1
    jHi = jLo+NBASE-1
  end if
else
#endif
  jLo = 1
  jHi = NIS
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_DISTRIBUTION
