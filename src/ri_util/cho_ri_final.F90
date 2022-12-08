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

subroutine Cho_RI_Final(irc,nVec_RI,l_nVec_RI)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_nVec_RI, nVec_RI(l_nVec_RI)
#include "cholesky.fh"

if (l_nVec_RI < nSym) then
  irc = 1
else
  irc = 0
  call Put_iArray('nVec_RI',nVec_RI,nSym)
end if

return

end subroutine Cho_RI_Final
