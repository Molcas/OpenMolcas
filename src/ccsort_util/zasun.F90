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

subroutine zasun(i1,length,valn,jn,kn,ln)
! control routine over zasun process
!
! i1 - number of pivot index (I)
! length - number of valid integrals in block (I)

use ccsort_global, only: mbas, nsize, zrkey
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: i1, length, jn(nsize,mbas), kn(nsize,mbas), ln(nsize,mbas)
real(kind=wp), intent(_IN_) :: valn(nsize,mbas)

if (zrkey == 1) then
  call zasun_zr(i1,length,valn,jn,kn,ln)
else
  call zasun_pck(i1,length,valn,jn,kn,ln)
end if

return

end subroutine zasun
