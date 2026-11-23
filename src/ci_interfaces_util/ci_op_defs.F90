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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The constants of the CI operations (CI_Operations), shared with the CI
!> libraries that it calls, and the CI library selected for the run (one
!> selection for every program: RASSCF and MCLR set it with CI_Lib_Select,
!> everyone else reads CI_Lib or asks CI_Lib_Is)
module CI_Op_Defs

use Definitions, only: iwp, u6

implicit none
private

!> The CI libraries (Lucia, SGUGA and Faroald with the values of RASSCF's
!> CI_Basis(); always compare against these names). CI_Lib_Lucia is
!> RASSCF's Lucia, CI_Lib_MCLR_Lucia the copy of Lucia in MCLR.
integer(kind=iwp), parameter :: CI_Lib_None = 0, CI_Lib_Lucia = 1, CI_Lib_SGUGA = 2, CI_Lib_Faroald = 3, CI_Lib_MCLR_Lucia = 4

!> Kinds of diagonal (CI_Diag)
integer(kind=iwp), parameter :: CI_Diag_Cnf_Avg = 1, CI_Diag_Exact = 2

!> Directions of the basis conversion: CI_To_Basis, stored basis -> basis
!> of the library; CI_To_File, the reverse
integer(kind=iwp), parameter :: CI_To_File = 0, CI_To_Basis = 1

!> The selected CI library (CI_Lib_*), written only by CI_Lib_Select. A
!> program that uses a CI library has to select it: the default is none.
!> Closing a library (CI_Close, CI_Lib_Close) frees its data but keeps the
!> selection.
integer(kind=iwp), protected :: CI_Lib = CI_Lib_None

public :: CI_Diag_Cnf_Avg, CI_Diag_Exact, CI_Lib, CI_Lib_Faroald, CI_Lib_Is, CI_Lib_Lucia, CI_Lib_MCLR_Lucia, CI_Lib_None, &
          CI_Lib_Select, CI_Lib_SGUGA, CI_To_Basis, CI_To_File

contains

!> Select the CI library of the run (CI_Lib_*; CI_Lib_None clears it).
!> Whether a library can do a given CI is decided by the program before
!> (RASSCF, proc_inp) or at the initialization of the library (MCLR,
!> CI_Lib_Initialize).
subroutine CI_Lib_Select(iLib)
  integer(kind=iwp), intent(in) :: iLib
  select case (iLib)
    case (CI_Lib_None,CI_Lib_Lucia,CI_Lib_SGUGA,CI_Lib_Faroald,CI_Lib_MCLR_Lucia)
      CI_Lib = iLib
    case default
      write(u6,'(A,I4)') ' CI_Lib_Select: unknown CI library',iLib
      call Abend()
  end select
end subroutine CI_Lib_Select

!> Is iLib the selected CI library?
function CI_Lib_Is(iLib) result(Is)
  integer(kind=iwp), intent(in) :: iLib
  logical(kind=iwp) :: Is
  Is = (CI_Lib == iLib)
end function CI_Lib_Is

end module CI_Op_Defs
