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

function MorsCre(IMORS,IPOS)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: MorsCre
integer(kind=iwp) :: IMORS, IPOS
integer(kind=iwp) :: ISGN
integer(kind=iwp), external :: MorsParity

MorsCre = 999999
if (btest(IMORS,IPOS-1)) return
MorsCre = ibset(IMORS,IPOS-1)
ISGN = MorsParity(IMORS/(2**(IPOS-1)))
MorsCre = ISGN*MorsCre

return

end function MorsCre
