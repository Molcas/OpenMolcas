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

subroutine mk_IOFF(IOFF,mSYM,NBASF,ISY12)

use Symmetry_Info, only: MUL
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mSYM, NBASF(mSym), ISY12
integer(kind=iwp), intent(out) :: IOFF(mSYM)
integer(kind=iwp) :: IOF, ISY1, ISY2, NB1, NB12, NB2

! FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF TDMSCR
IOF = 0
IOFF(:) = 0
do ISY1=1,mSYM
  ISY2 = MUL(ISY1,ISY12)
  if (ISY1 < ISY2) cycle
  IOFF(ISY1) = IOF
  IOFF(ISY2) = IOF
  NB1 = NBASF(ISY1)
  NB2 = NBASF(ISY2)
  NB12 = NB1*NB2
  if (ISY1 == ISY2) NB12 = (NB12+NB1)/2
  IOF = IOF+NB12
end do

end subroutine mk_IOFF
