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

subroutine GT1DIA(H1DIA)
! Obtain diagonal of one electron matrix over active orbitals
!
! Dec 97 : obtained from INT1O

use lucia_data, only: IBSO, INT1O, IREOTS, ISMFTO, NACOB, PINT1
use Definitions, only: wp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: H1DIA(*)

! Global pointers

!INA call GT1DIS(H1DIA,IREOTS(1+NINOB),PINT1,INT1,ISMFTO,IBSO,NACOB)
call GT1DIS(H1DIA,IREOTS,PINT1,INT1O,ISMFTO,IBSO,NACOB)

end subroutine GT1DIA
