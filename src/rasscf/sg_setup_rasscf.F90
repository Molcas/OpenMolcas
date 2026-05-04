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

#ifdef _DMRG_
subroutine SG_Setup_RASSCF(DBG,SkipGUGA,initial_occ)
#else
subroutine SG_Setup_RASSCF(DBG,SkipGUGA)
#endif

use fciqmc, only: DoNECI
use fcidump, only: DumpOnly
use CC_CI_mod, only: Do_CC_CI
use gas_data, only: iDoGAS
use rasscf_global, only: DoBlockDMRG
use general_data, only: nSym, nActel, iSpin, nHole1, nElec3, nRs1, nRs2, nRs3, STSYM, nCOnf
#ifdef _DMRG_
use rasscf_global, only: DoDMRG
use input_ras, only: Key
use stdalloc, only: mma_deallocate
#endif
use sguga, only: CIS, EXS, SGS, MkISM_RASSCF
use definitions, only: u6, iwp, wp

#ifdef _DMRG_
integer(kind=iwp), allocatable, intent(inout) :: initial_occ(:,:)
#endif

logical(kind=iwp), intent(inout):: DBG,SkipGUGA
real(kind=wp) Eterna_1, Eterna_2, dum1, dum2, dum3

! Construct the Guga tables

if (.not. (DoNECI .or. Do_CC_CI .or. DumpOnly .or. SkipGUGA)) then
  ! right now skip most part of gugactl for GAS, but only call mkism.
  if (.not. iDoGas) then
    ! DMRG calculation no need the SG_Init_RASSCF subroutine
#   ifdef _DMRG_
    if (Key('DMRG') .or. doDMRG) then
      call mma_deallocate(initial_occ)
      SkipGUGA = .true.
    else
#   endif
      call Timing(Eterna_1,dum1,dum2,dum3)
      if (DBG) write(u6,*) ' Call GugaCtl'
      call SG_Init_RASSCF(nSym,nActEl,iSpin,               &
                          SGS,CIS,EXS,                     &
                          nHole1,nElec3,nRs1,nRs2,nRs3,    &
                          STSYM,DoBlockDMRG)
!     (SGS%IFRAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
!     IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
!     SET UP THE REINDEXING TABLE
      call SETSXCI()
      NCONF = CIS%NCSF(STSYM)

      call Timing(Eterna_2,dum1,dum2,dum3)
#   ifdef _DMRG_
    end if
#   endif
  else  ! if iDoGas
    call mkism_rasscf(SGS)
  end if
end if

end subroutine sg_setup_rasscf
