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

subroutine inisewm(prgnam,ndiff)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: prgnam
integer(kind=iwp), intent(inout) :: nDiff
character(len=16) :: pgnm_local
logical(kind=iwp) :: DoRys

pgnm_local = prgnam
call locase(pgnm_local)

select case (pgnm_local)
  case ('seward','slapaf')
  case ('scf','mltpl','alaska','mckinley','espf')
    DoRys = .true.  ! for Schwarz prescreening
    call inisew(DoRys,ndiff)
  case default
    DoRys = .false.
    call inisew(DoRys,ndiff)
end select

return

end subroutine inisewm
