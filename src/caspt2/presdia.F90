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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
!> @brief Apply the resolvent of the diagonal of H0-E0 to a vector
!> @author Per-&Aring;ke Malmqvist, 1994
!>
!> @details
!> ::presdia stands for "parallel resolvent diagonal" and it applies
!> the resolvent of the diagonal of \f$ (H_0 - E_0) \f$ to a coefficient
!> vector stored at position \p IVEC on \p LUSOLV. The result is stored
!> in the vector at position \p JVEC. Furthermore, it computes the
!> "overlap" as \p IVEC squared divided by the resolvent.
!> The parallel part is taken care of in subroutine calls \p RHS_.
!> Potential shifts and modification of \f$ H_0 \f$ are considered
!> in ::rhs_resdia
!>
!> @param[in]    IVEC   Vector position to which the res is applied
!> @param[in]    JVEC   Vector position where the result is saved
!> @param[out]   OVLAPS Array containing the overlaps

subroutine PRESDIA(IVEC,JVEC,OVLAPS)

use EQSOLV, only: IDBMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: MxCase, nASup, nInDep, nISup, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC
real(kind=wp), intent(inout) :: OVLAPS(0:8,0:MXCASE)
integer(kind=iwp) :: ICASE, ISYM, JD, lg_V, NAS, NIN, NIS
real(kind=wp) :: DOVL, OVL, OVLSUM, OVLTOT
real(kind=wp), allocatable :: BD(:), ID(:)

! Apply the resolvent of the diagonal part of H0 to a coefficient
! vector in vector nr. IVEC on LUSOLV. Put the results in vector
! nr. JVEC. Also compute overlaps, see OVLVEC for structure.

OVLTOT = Zero
OVLAPS(:,:) = Zero

do ICASE=1,13
  OVLSUM = Zero
  do ISYM=1,NSYM
    OVL = Zero
    NIN = NINDEP(ISYM,ICASE)
    if (NIN <= 0) cycle
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    ! Remember: NIN values in BDIAG, but must read NAS for correct
    ! positioning.
    call mma_allocate(BD,NAS,LABEL='BD')
    call mma_allocate(ID,NIS,LABEL='ID')

    ! Read the B matrix of each case
    JD = IDBMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,2,BD,NAS,JD)
    call DDAFILE(LUSBT,2,ID,NIS,JD)

    call RHS_ALLO(NIN,NIS,lg_V)
    call RHS_READ(NIN,NIS,lg_V,ICASE,ISYM,IVEC)
    call RHS_RESDIA(NIN,NIS,lg_V,BD,ID,DOVL)

    OVL = OVL+DOVL

    call RHS_SAVE(NIN,NIS,lg_V,ICASE,ISYM,JVEC)
    call RHS_FREE(lg_V)

    call mma_deallocate(BD)
    call mma_deallocate(ID)
    OVLAPS(ISYM,0) = OVLAPS(ISYM,0)+OVL
    OVLSUM = OVLSUM+OVL
  end do
  OVLAPS(0,ICASE) = OVLSUM
  OVLTOT = OVLTOT+OVLSUM
end do
OVLAPS(0,0) = OVLTOT

end subroutine PRESDIA
