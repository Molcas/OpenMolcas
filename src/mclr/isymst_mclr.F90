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

integer function ISYMST_MCLR(STRING,NEL)
! Master routine for symmetry of string

use MCLR_Data, only: ISMFTO

implicit none
! Specific input
integer STRING(*), NEL
integer ISYM, IEL, JSYM, IVV, KVV

ISYM = 1
do IEL=1,NEL
  JSYM = ISMFTO(STRING(IEL))-1
  IVV = ISYM-1
  KVV = ieor(IVV,JSYM)
  ISYM = KVV+1
end do
ISYMST_MCLR = ISYM

end function ISYMST_MCLR
