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

function MorsSymm(IMORS,ISMARR)

use Symmetry_Info, only: MUL
use Cntrl, only: MORSBITS
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) MorsSymm
integer(kind=iwp) IMORS, ISMARR(*)
integer(kind=iwp) I, IB, IBIT

MorsSymm = 1
if (IMORS < 0) then
  write(u6,*) ' MorsSymm: Bad IMORS=',IMORS
  call ABEND()
end if
IB = IMORS
do I=1,MORSBITS
  IBIT = mod(IB,2)
  IB = IB/2
  if (IBIT == 1) MorsSymm = MUL(MorsSymm,ISMARR(I))
end do

end function MorsSymm
