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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA(rc)

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Driver routine for the calculation of correlation energies in the
! random-phase approximation (RPA) using Cholesky/DF integrals.
!
! NOTE: conventional integrals are not implemented!

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: rc
#include "warnings.h"
integer(kind=iwp) :: irc
character(len=80) :: string
character(len=*), parameter :: SecNam = 'RPA'

!=======================================================================
! Enter
!=======================================================================

rc = _RC_ALL_IS_WELL_ ! init return code
irc = 0 ! init internal return code

!=======================================================================
! Setup: initialize data, read reference orbitals and process input
!=======================================================================

call StatusLine('RPA: ','Setup')
call RPA_Setup()

!=======================================================================
! Exit after cleanup
!=======================================================================

call StatusLine('RPA: ','Cleanup')
call RPA_Cleanup(irc)
if (irc /= 0) then
  write(string,'(A,A,I4)') SecNam,': Cleanup failed! rc=',irc
  call WarningMessage(2,string)
  if (rc == _RC_ALL_IS_WELL_) rc = _RC_INTERNAL_ERROR_
end if

end subroutine RPA
