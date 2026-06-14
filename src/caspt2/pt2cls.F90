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

subroutine PT2CLS()

use INPUTDATA, only: CLEANUP_INPUT
use SUPERINDEX, only: SUPFREE
use PT2WFN, only: PT2WFN_CLOSE
use sguga, only: SG_Free
#ifdef _DMRG_
use qcmaquis_interface, only: qcmaquis_interface_deinit
use qcmaquis_interface_cfg, only: dmrg_file
use qcmaquis_info, only: qcmaquis_info_deinit
use caspt2_module, only: DMRG
#endif
#if 0
! NOT TESTED
use OFembed, only: FMaux
#endif
use ChoCASPT2, only: NASplit, NISplit, NumCho_PT2
use caspt2_global, only: CMOPT2, DMIX, DREF, DWGT, FIFA, FIMO, IDCIEX, IDSCT, IDTCEX, PREF, TAT, TORB, Weight, SGS, CIS, TRS, EXS
use caspt2_module, only: IfChol, nAsh, nIsh, nSsh, nSym
use stdalloc, only: mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: irc, iSym

if (IfChol) then
  ! Finalize Cholesky information
  call Cho_X_Final(irc)
  if (irc /= 0) then
    write(u6,*) 'CASPT2: Non-zero rc in Cho_X_Final'
    call QUIT(irc)
  end if
  ! Delete files of MO Cholesky vectors.
  do iSym=1,nSym
    call Cho_Caspt2_OpenF(3,1,iSym,nIsplit(iSym)) !close+delete
    call Cho_Caspt2_OpenF(3,2,iSym,nAsplit(iSym)) !close+delete
  end do
  ! Deallocate memory used as index arrays in the setup section
  call setup_cho(nSym,nIsh,nAsh,nSsh,NumCho_pt2,'FREE')
# if 0
  ! NOT TESTED
  call mma_deallocate(FMaux,safe='*')
# endif
  ! deallocate chovec_io arrays
  call trachosz_free()
end if

! Deallocate SGUGA tables:
call SG_Free(SGS,CIS,TRS,EXS)

! dealloacte DMRG stuff
#ifdef _DMRG_
if (DMRG) then
  call mma_deallocate(dmrg_file%qcmaquis_checkpoint_file)
  call qcmaquis_info_deinit()
  call qcmaquis_interface_deinit()
end if
#endif

! Deallocate MAGEB, etc, superindex tables:
call SUPFREE()
! Deallocate global array for Fock matrix, etc:
call mma_deallocate(IDCIEX)
call mma_deallocate(IDTCEX)
call mma_deallocate(FIFA)
call mma_deallocate(FIMO)
call mma_deallocate(DREF)
call mma_deallocate(PREF)
call mma_deallocate(DMIX)
call mma_deallocate(DWGT)
! Deallocate global orbital transformation arrays:
call mma_deallocate(TORB)
call mma_deallocate(TAT)
! Deallocate global orbital arrays:
call mma_deallocate(CMOPT2)
! Deallocate global RHS disk offsets (allocated in eqctl1):
call mma_deallocate(IDSCT)

call pt2wfn_close()
! Close all files:
call CLSFLS_CASPT2()

! Weight array is allocated in refwfn_info
call mma_deallocate(Weight)

! free input struct
call CleanUp_Input()

end subroutine PT2CLS
