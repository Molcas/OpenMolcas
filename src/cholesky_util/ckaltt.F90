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

subroutine CkAltT(K_Lap,R,T,NG)
!-----------------------------------------------------------------------
! Function : Check the alternation of array T
!-----------------------------------------------------------------------

use ReMez_mod

implicit real*8(A-H,O-Z)
parameter(ONE=1.0D+00)
real*8 T(40)
logical NG

IDimEnd = 2*K_Lap+1
NG = .false.

do I=1,IDimEnd
  if (I == 1) then
    A = ONE
  else
    A = T(I-1)
  end if
  if (I == IDimEnd) then
    B = R
  else
    B = T(I)
  end if
  if (A >= B) then
    write(IW,'(A,I3)') 'The sign of T is wrong at I =',I
    NG = .true.
    goto 999
  end if
end do

999 continue
return

end subroutine CkAltT
