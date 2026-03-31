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
integer(kind=iwp) :: ISGN, MASK
integer(kind=iwp), external :: MorsParity

MorsAnn = 999999
MASK = 2**(IPOS-1)
if (iand(MASK,IMORS) == 0) return
MorsAnn = IMORS-MASK
ISGN = MorsParity(MorsAnn/MASK)
MorsAnn = ISGN*MorsAnn

return

end function MorsAnn
