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

!> The constants of the CI operations (CI_Operations), shared with the CI
!> libraries that it calls
module CI_Op_Defs

use Definitions, only: iwp

implicit none
private

!> The CI libraries (Lucia, SGUGA and Faroald with the values of RASSCF's
!> CI_Basis(); always compare against these names)
integer(kind=iwp), parameter :: CI_Lib_None = 0, CI_Lib_Lucia = 1, CI_Lib_SGUGA = 2, CI_Lib_Faroald = 3, CI_Lib_MCLR_Lucia = 4

!> Kinds of diagonal (CI_Diag)
integer(kind=iwp), parameter :: CI_Diag_Cnf_Avg = 1, CI_Diag_Exact = 2

!> Directions of the basis conversion: CI_To_Basis, stored basis -> basis
!> of the library; CI_To_File, the reverse
integer(kind=iwp), parameter :: CI_To_File = 0, CI_To_Basis = 1

public :: CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib_Faroald, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_None, CI_Lib_SGUGA, CI_To_Basis, &
          CI_To_File

end module CI_Op_Defs
