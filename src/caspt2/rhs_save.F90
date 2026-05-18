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

subroutine RHS_SAVE(NIN,NIS,lg_W,iCASE,iSYM,iVEC)
!SVC: this routine reads an RHS array in SR format from disk

use fake_GA, only: GA_Arrays
use caspt2_global, only: LURHS
use caspt2_module, only: IOFFRHS
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use Definitions, only: u6
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NIN, NIS, lg_W, iCASE, iSYM, iVEC
integer(kind=iwp) :: IDISK, NW
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: IEND, ISTA, JEND, JSTA, LDW, mpt_W, myRank, NWPROC
#include "global.fh"
#include "mafdecls.fh"

if (Is_Real_Par()) then
  call GA_Sync()
  myRank = GA_NodeID()
  call GA_Distribution(lg_W,myRank,ISTA,IEND,JSTA,JEND)
  if ((IEND-ISTA+1 == NIN) .and. (ISTA > 0)) then
    call GA_Access(lg_W,ISTA,IEND,JSTA,JEND,mpt_W,LDW)
    if (LDW /= NIN) then
      write(u6,*) 'RHS_SAVE: Assumption NIN==LDW wrong'
      call AbEnd()
    end if
    NWPROC = NIN*(JEND-JSTA+1)
    IDISK = IOFFRHS(ISYM,ICASE)
    call DDAFILE(LURHS(IVEC),1,DBL_MB(mpt_W),NWPROC,IDISK)
    call GA_Release(lg_W,ISTA,IEND,JSTA,JEND)
  end if
  call GA_Sync()
else
#endif
  NW = NIN*NIS
  IDISK = IOFFRHS(ISYM,ICASE)
  call DDAFILE(LURHS(IVEC),1,GA_Arrays(lg_W)%A,NW,IDISK)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_SAVE
