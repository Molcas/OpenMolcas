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

subroutine RHS_READ(NIN,NIS,lg_W,iCASE,iSYM,iVEC)
!SVC: this routine reads an RHS array in SR format from disk

use definitions, only: iwp
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use caspt2_global, only: LURHS
use fake_GA, only: GA_Arrays
use caspt2_module, only: IOFFRHS

implicit none
integer(kind=iwp), intent(in) :: NIN, NIS, lg_W, iCASE, iSYM, iVEC
integer(kind=iwp) IDISK, NWPROC
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
integer(kind=iwp) myRank, ISTA, IEND, JSTA, JEND, mpt_W, LDW

if (Is_Real_Par()) then
  call GA_Sync()
  myRank = GA_NodeID()
  call GA_Distribution(lg_W,myRank,ISTA,IEND,JSTA,JEND)
  if ((IEND-ISTA+1 == NIN) .and. (ISTA > 0)) then
    call GA_Access(lg_W,ISTA,IEND,JSTA,JEND,mpt_W,LDW)
    if (LDW /= NIN) then
      write(6,*) 'RHS_READ: Assumption NIN==LDW wrong'
      call AbEnd()
    end if
    NWPROC = NIN*(JEND-JSTA+1)
    IDISK = IOFFRHS(ISYM,ICASE)
    call DDAFILE(LURHS(IVEC),2,DBL_MB(mpt_W),NWPROC,IDISK)
    call GA_Release_Update(lg_W,ISTA,IEND,JSTA,JEND)
  end if
  call GA_Sync()
else
#endif
  NWPROC = NIN*NIS
  IDISK = IOFFRHS(ISYM,ICASE)
  call DDAFILE(LURHS(IVEC),2,GA_Arrays(lg_W)%A,NWPROC,IDISK)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_READ
