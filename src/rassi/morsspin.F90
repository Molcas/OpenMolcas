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

function MorsSpin(IMORS,MS2ARR)

use Cntrl, only: MORSBITS
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: MorsSpin
integer(kind=iwp) :: IMORS, MS2ARR(*)
integer(kind=iwp) :: I, IB, IBIT

MorsSpin = 0
if (IMORS < 0) then
  write(u6,*) ' MorsSpin: Bad IMORS=',IMORS
  call ABEND()
end if
IB = IMORS
do I=1,MORSBITS
  IBIT = mod(IB,2)
  IB = IB/2
  if (IBIT == 1) MorsSpin = MorsSpin+MS2ARR(I)
end do

end function MorsSpin
