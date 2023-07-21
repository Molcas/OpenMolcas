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

use ReMez_mod

implicit real*8(A-H,O-Z)
parameter(ZERO=0.0D+00)
real*8 DD(82)
logical NG

IDimEnd = 2*K_Lap
NG = .false.

do I=1,IDimEnd
  A = DD(I)*DD(I+1)
  if (A >= ZERO) then
    NG = .true.
    write(IW,'(A,I3)') 'DD sign is wrong at I =',I
  end if
end do

return

end subroutine CkAltD
