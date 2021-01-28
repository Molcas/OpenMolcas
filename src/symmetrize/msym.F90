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
! Copyright (C) 2015, Marcus Johansson                                 *
!***********************************************************************

subroutine fmsym(rc)

use Definitions, only: wp, iwp
implicit none
integer(kind=iwp), intent(out) :: rc
real(kind=wp) :: ctx

call fmsym_create_context(ctx)
call fmsym_set_elements(ctx)
call fmsym_find_symmetry(ctx)
call fmsym_symmetrize_molecule(ctx)
call fmsym_generate_orbital_subspaces(ctx)
call fmsym_symmetrize_orb_file(ctx,'INPORB')
call fmsym_release_context(ctx)
rc = 0

return

end subroutine fmsym
