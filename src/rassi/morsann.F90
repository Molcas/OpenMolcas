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

function MorsAnn(IMORS,IPOS)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: MorsAnn
integer(kind=iwp) :: IMORS, IPOS
integer(kind=iwp) :: ISGN
integer(kind=iwp), external :: MorsParity

MorsAnn = 999999
if (.not. btest(IMORS,IPOS-1)) return
MorsAnn = ibclr(IMORS,IPOS-1)
ISGN = MorsParity(MorsAnn/(2**(IPOS-1)))
MorsAnn = ISGN*MorsAnn

return

end function MorsAnn
