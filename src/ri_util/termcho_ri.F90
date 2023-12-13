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

subroutine TermCho_RI(irc,nVec_RI,l_nVec_RI)

use Cholesky, only: MySP
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_nVec_RI, nVec_RI(l_nVec_RI) ! #RI vectors per irrep on this node

irc = 0

! Save number of vectors and other info on runfile.
! -------------------------------------------------

call Cho_Final(.false.)
call Cho_RI_Final(irc,nVec_RI,l_nVec_RI)
if (irc /= 0) return

! Close storage files.
! --------------------

call Cho_P_OpenVR(2)

! Deallocate index arrays.
! ------------------------

call Cho_X_Dealloc(irc)
if (irc /= 0) return

! More deallocations.
! -------------------

if (allocated(MySP)) call mma_deallocate(MySP)

end subroutine TermCho_RI
