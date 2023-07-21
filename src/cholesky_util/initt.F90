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

subroutine InitT(K_Lap,T,RIni,RNew)
!-----------------------------------------------------------------------
! Function : Set the initial T()
!-----------------------------------------------------------------------

use ReMez_mod

implicit real*8(A-H,O-Z)
parameter(ONE=1.0D+00)
real*8 T(40)
logical Dbg

Dbg = .false.
Val = (RNew-ONE)/(RIni-ONE)

if (Dbg) write(IW,'(A)') 'T(I)'
do I=1,2*K_Lap
  T(I) = ONE+(T(I)-ONE)*Val
  if (Dbg) write(IW,'(I3,2X,F20.14)') I,T(I)
end do

return

end subroutine InitT
