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

function ISYMST_MCLR(STRING,NEL)
! Master routine for symmetry of string

use Symmetry_Info, only: Mul
use MCLR_Data, only: ISMFTO
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISYMST_MCLR
integer(kind=iwp), intent(in) :: NEL, STRING(NEL)
integer(kind=iwp) :: IEL

ISYMST_MCLR = 1
do IEL=1,NEL
  ISYMST_MCLR = Mul(ISYMST_MCLR,ISMFTO(STRING(IEL)))
end do

end function ISYMST_MCLR
