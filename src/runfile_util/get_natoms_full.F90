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
! Copyright (C) 2015, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Get_nAtoms_Full
!
!> @author I. Fdez. Galv&aacute;n
!>
!> @details
!> Get number of all atoms (not only symmetry unique) from RUNFILE,
!> including MM atoms otherwise invisible to gateway/slapaf.
!>
!> @param[out] nAtoms_Full Number of all atoms in the system
!***********************************************************************

subroutine Get_nAtoms_Full(nAtoms_Full)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: nAtoms_Full
integer(kind=iwp) :: nAtMM, nAtom
logical(kind=iwp) :: Found

call Get_nAtoms_All(nAtom)
call Qpg_dArray('MMO Coords',Found,nAtMM)
nAtoms_Full = nAtom+nAtMM/3

return

end subroutine Get_nAtoms_Full
