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

subroutine PSGMDIA(ALPHA,BETA,IVEC,JVEC)
! Compute |JVEC> := BETA*|JVEC> + ALPHA*(H0(diag)-E0)*|IVEC>
! If real_shift /= Zero or imag_shift /= Zero, use a modified H0

use definitions, only: iwp, wp
use constants, only: Zero
use caspt2_global, only: LUSBT
use EQSOLV, only: IDBMat
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: nSym, nInDep, nASup, nISup

implicit none
real(kind=wp), intent(in) :: ALPHA, BETA
integer(kind=iwp), intent(in) :: IVEC, JVEC
integer(kind=iwp) ICASE, ISYM, NIN, NAS, NIS, JD, lg_V1, lg_V2
real(kind=wp), allocatable :: BD(:), ID(:)

do ICASE=1,13
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    if (NIN == 0) cycle
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NIS == 0) cycle
    ! Remember: NIN values in BDIAG, but must read NAS for correct
    ! positioning.
    call mma_allocate(BD,NAS,LABEL='BD')
    call mma_allocate(ID,NIS,LABEL='ID')
    JD = IDBMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,2,BD,NAS,JD)
    call DDAFILE(LUSBT,2,ID,NIS,JD)

    call RHS_ALLO(NIN,NIS,lg_V2)

    if (BETA /= Zero) then
      call RHS_READ(NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
      if (BETA /= 1) then
        call RHS_SCAL(NIN,NIS,lg_V2,BETA)
      !else
      !  call RHS_SCAL(NIN,NIS,lg_V2,Zero)
      end if
    end if

    if (ALPHA /= Zero) then
      if (BETA /= Zero) then
        call RHS_ALLO(NIN,NIS,lg_V1)
        call RHS_READ(NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
        call RHS_SGMDIA(NIN,NIS,lg_V1,BD,ID)
        call RHS_DAXPY(NIN,NIS,ALPHA,lg_V1,lg_V2)
        call RHS_FREE(lg_V1)
      else
        call RHS_READ(NIN,NIS,lg_V2,ICASE,ISYM,IVEC)
        call RHS_SGMDIA(NIN,NIS,lg_V2,BD,ID)
        call RHS_SCAL(NIN,NIS,lg_V2,ALPHA)
      end if
    end if

    call RHS_SAVE(NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
    call RHS_FREE(lg_V2)
    call mma_deallocate(BD)
    call mma_deallocate(ID)
  end do
end do

end subroutine PSGMDIA
