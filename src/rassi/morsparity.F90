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

function MorsParity(IMORS)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: MorsParity
integer(kind=iwp) :: IMORS
integer(kind=iwp) :: I1, I2, I3, J
integer(kind=iwp), parameter :: ISG(0:15) = [1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1]

MorsParity = 0 ! dummy initialize
if (IMORS >= 0) then
  I1 = IMORS/16
  J = IMORS-16*I1
  MorsParity = ISG(J)
  if (I1 == 0) return
  I2 = I1/16
  J = I1-16*I2
  MorsParity = MorsParity*ISG(J)
  if (I2 == 0) return
  I3 = I2/16
  J = I2-16*I3
  MorsParity = MorsParity*ISG(J)
  if (I3 == 0) return
end if
write(u6,*) ' MorsParity: Bad IMORS=',IMORS
call ABEND()

end function MorsParity
