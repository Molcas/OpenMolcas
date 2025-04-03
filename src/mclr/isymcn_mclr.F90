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

integer function ISYMCN_MCLR(IOP,NOPEN)
! Master routine for symmetry of configuration
! with NOPEN singly occupied shells

use Symmetry_Info, only: Mul
use MCLR_Data, only: ISMFTO

implicit none
integer IOP(*)
integer NOPEN
integer IEL

ISYMCN_MCLR = 1
do IEL=1,NOPEN
  ISYMCN_MCLR = Mul(ISYMCN_MCLR,ISMFTO(IOP(IEL)))
end do

end function ISYMCN_MCLR
