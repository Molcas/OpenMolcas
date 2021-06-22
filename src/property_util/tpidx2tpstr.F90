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

subroutine tpidx2tpstr(TYPEINDEX,TYPESTRING,NB)

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: TYPEINDEX(*), NB
character, intent(_OUT_) :: TYPESTRING(*)
integer(kind=iwp) :: i

do i=1,NB
  select case (TYPEINDEX(i))
    case (1)
      TYPESTRING(i) = 'F'
    case (2)
      TYPESTRING(i) = 'I'
    case (3)
      TYPESTRING(i) = '1'
    case (4)
      TYPESTRING(i) = '2'
    case (5)
      TYPESTRING(i) = '3'
    case (6)
      TYPESTRING(i) = 'S'
    case (7)
      TYPESTRING(i) = 'D'
  end select
end do

end subroutine tpidx2tpstr
