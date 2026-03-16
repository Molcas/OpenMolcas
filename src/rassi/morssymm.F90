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

integer function MorsSymm(IMORS,ISMARR)

use Symmetry_Info, only: MUL
use Definitions, only: u6

implicit none
integer MORSBITS
parameter(MORSBITS=8)
dimension ISMARR(*)
integer I, IB, IBIT, IMORS, ISMARR

MorsSymm = 1
if (IMORS < 0) goto 99
IB = IMORS
do I=1,MORSBITS
  IBIT = mod(IB,2)
  IB = IB/2
  if (IBIT == 1) MorsSymm = MUL(MorsSymm,ISMARR(I))
end do

return
99 continue
write(u6,*) ' MorsSymm: Bad IMORS=',IMORS
call ABEND()

end function MorsSymm
