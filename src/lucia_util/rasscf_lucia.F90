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

module RASSCF_LUCIA

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: ini_h0, kvec3_length = 0, Memory_Needed_Lucia = 0
logical(kind=iwp) :: Sigma_on_disk = .false.
real(kind=wp), allocatable :: DStmp(:), Dtmp(:), PAtmp(:), Pscr(:), Ptmp(:), RF1(:), RF2(:)

public :: ini_h0, kvec3_length, Memory_Needed_Lucia, Sigma_on_disk, DStmp, Dtmp, PAtmp, Pscr, Ptmp, RF1, RF2

end module RASSCF_LUCIA
