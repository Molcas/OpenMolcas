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
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine active_classical_coulomb(DmatDmat, eri, coulomb_energy)
  use stdalloc, only: mma_allocate, mma_deallocate
  use definitions, only: wp
  implicit none

  real(kind=wp), dimension(*), intent(in) :: DmatDmat, eri
  real(kind=wp), intent(out) :: coulomb_energy

#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"

  real(kind=wp), allocatable, dimension(:) :: q
  integer :: iSym, nActive, NUVX, nOrbitals, nInactive
  integer :: ISTP, jstf, nt, ntt

  coulomb_energy = 0.0D0

  do iSym=1, nSym
    nActive = nAsh(iSym)
    nOrbitals = nOrb(iSym)
    nInactive = NISH(ISYM)

    if (nActive /= 0) then
      ISTP=ISTORP(ISYM)+1
      JSTF=ISTORD(ISYM)+1
      NUVX = (ISTORP(ISYM+1)-ISTORP(ISYM))/nActive

      call mma_allocate(q, nOrbitals*nActive, "q")
      CALL DGEMM_('N','N', nOrbitals,nActive, NUVX, 1.0d0, eri(JSTF), &
              nOrbitals, DmatDmat(ISTP), NUVX, 0.0d0, Q, nOrbitals)

      do nt=1, nActive
        ntt = (nt-1)*nOrbitals + nInactive + nt
        coulomb_energy = coulomb_energy + q(ntt)
      end do

      call mma_deallocate(q)
    end if
  end do

  coulomb_energy = 0.5D0 * coulomb_energy
end subroutine active_classical_coulomb