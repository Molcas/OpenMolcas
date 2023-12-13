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

subroutine tpstr2tpidx(TYPESTRING,TYPEINDEX,NB)

use Definitions, only: iwp

#include "intent.fh"

implicit none
character, intent(in) :: TYPESTRING(*)
integer(kind=iwp), intent(_OUT_) :: TYPEINDEX(*)
integer(kind=iwp), intent(in) :: NB
integer(kind=iwp) :: i

do i=1,NB
  select case (TYPESTRING(i))
    case ('F','f')
      TYPEINDEX(i) = 1
    case ('I','i')
      TYPEINDEX(i) = 2
    case ('1')
      TYPEINDEX(i) = 3
    case ('2')
      TYPEINDEX(i) = 4
    case ('3')
      TYPEINDEX(i) = 5
    case ('S','s')
      TYPEINDEX(i) = 6
    case ('D','d')
      TYPEINDEX(i) = 7
  end select
end do

end subroutine tpstr2tpidx
