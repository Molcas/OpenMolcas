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

subroutine get_ref_energy(nroots, jobold, iadr19, ref_energy)
  ! Get reference energies from JOBOLD and stores them into ref_energy
  use definitions, only: wp
  implicit none

  integer, intent(in) :: nroots, jobold
  integer, dimension(15), intent(in) :: iadr19
  real(kind=wp), dimension(nroots), intent(out) :: ref_energy
! for mxiter, mxroot
#include "rasdim.fh"

  integer :: jdisk, iter, root, maybe
  real(kind=wp) :: aemax
  real(kind=wp), dimension(mxiter*mxroot) :: elist

  jdisk = iAdr19(6)
  call DDaFile(JOBOLD, 2, elist, MXROOT*MXITER, jdisk)

  maybe = 0
  do iter=1, mxiter
    aemax = 0.0D0
    do root=1, mxroot
      AEMAX=max(aemax, abs(elist(mxroot*(iter-1) + root)))
    end do
    if (abs(AEMAX) <= 1.0D-12) exit
    maybe=iter
  end do

  do root=1, nroots
    ref_energy(root) = elist(mxroot*(maybe-1)+root)
  end do

end subroutine get_ref_energy