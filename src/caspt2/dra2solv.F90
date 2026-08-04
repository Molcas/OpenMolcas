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

subroutine DRA2SOLV(NAS,NIS,iCASE,iSYM,iVEC)
!SVC: FIXME: this temporary routine copies the RHS arrays from DRAs to
!     LUSOLV and should be removed once the full parallelization is in
!     place and transition is no longer needed.

use caspt2_global, only: IDSCT, LUSOLV
use EQSOLV, only: MXSCT
use fake_GA, only: GA_Arrays
use caspt2_module, only: MXCASE
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Definitions, only: wp, u6
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NAS, NIS, iCASE, iSYM, iVEC
integer(kind=iwp) :: IDISK, lg_W
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: IEND, iMax, ISTA, NCOL, NW
real(kind=wp), allocatable :: TMPW(:)
#endif

!SVC: Read the global array from disk
call RHS_ALLO(NAS,NIS,lg_W)
call RHS_READ(NAS,NIS,lg_W,ICASE,ISYM,IVEC)

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  !SVC: Only the master process writes to LUSOLV!!
  !SVC: be careful to only call one-sided operations
  if (KING()) then
    !SVC: write the global array to LUSOLV
    !     when all routines have been adapted, this can go into SYNRHS and
    !     later it should be completely removed when VCUTIL and SGM are
    !     adapted for handling global + disk resident arrays then, also
    !     remove iCASE,iSYM,iVEC from the call, as they are no longer needed
    !     in that case.
    call mma_MaxDBLE(iMax)
    !-SVC: GA_Get does not like large buffer sizes, put upper limit at 1GB
    iMax = min(nint(0.95_wp*iMax),134217728)
    NCOL = min(iMAX,NAS*NIS)/NAS
    if (NCOL <= 0) then
      write(u6,*) 'Not enough memory in DRA2SOLV, aborting...'
      call AbEnd()
    end if
    NW = NAS*NCOL
    call mma_allocate(TMPW,NW,LABEL='TMPW')
    !SVC: Write local array to LUSOLV
    IDISK = IDSCT(1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
    do ISTA=1,NIS,NCOL
      IEND = min(ISTA+NCOL-1,NIS)
      call GA_Get(lg_W,1,NAS,ISTA,IEND,TMPW,NAS)
      call DDAFILE(LUSOLV,1,TMPW,NAS*(IEND-ISTA+1),IDISK)
    end do
    call mma_deallocate(TMPW)
  end if
  call GASync()
  !SVC: Destroy the global array
  !       bStat=GA_Destroy(lg_W)
else
#endif
  IDISK = IDSCT(1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
  call DDAFILE(LUSOLV,1,GA_Arrays(lg_W)%A,NAS*NIS,IDISK)
#ifdef _MOLCAS_MPP_
end if
#endif

call RHS_FREE(lg_W)

end subroutine DRA2SOLV
