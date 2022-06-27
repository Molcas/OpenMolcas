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

module loprop_arrays

use Definitions, only: wp, iwp

implicit none
private

! Encapsulating all these allocatable arrays in a derived type allows
! passing the type as an argument without an explicit interface, and
! allocating in the inner subroutine. It also allows having several
! simultaneous contexts with independent arrays (e.g. in a recursive
! calling sequence)

type LP_context_type
  integer(kind=iwp), allocatable :: ANr(:), center(:), otype(:)
  real(kind=wp), allocatable :: C(:,:), P(:,:), PInv(:,:), Q_Nuc(:)
end type LP_context_type

public :: LP_context_type

end module loprop_arrays
