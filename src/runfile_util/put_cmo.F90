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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************
!  Put_CMO
!
!> @brief
!>   Write the symmetry blocked MO coefficients on the runfile
!> @author R. Lindh
!>
!> @details
!> The utility will write the symmetry blocked MO coefficients on the runfile.
!>
!> @param[in] CMO  Array of symmetry blocked MO coefficients
!> @param[in] nCMO Number of elements in \p CMO
!***********************************************************************

subroutine Put_CMO(CMO,nCMO)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nCMO
real(kind=wp) :: CMO(nCMO)
character(len=24) :: Label

Label = 'Last orbitals'
call Put_dArray(Label,CMO,nCMO)

return

end subroutine Put_CMO
