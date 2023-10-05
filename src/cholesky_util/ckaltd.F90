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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine CkAltD(K_Lap,DD,NG)
!-----------------------------------------------------------------------
! Function : Check the alternation of array DD
!-----------------------------------------------------------------------

use ReMez_mod, only: IW
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: DD(82)
logical(kind=iwp), intent(out) :: NG
integer(kind=iwp) :: I, IDimEnd
real(kind=wp) :: A

IDimEnd = 2*K_Lap
NG = .false.

do I=1,IDimEnd
  A = DD(I)*DD(I+1)
  if (A >= Zero) then
    NG = .true.
    write(IW,'(A,I3)') 'DD sign is wrong at I =',I
  end if
end do

return

end subroutine CkAltD
