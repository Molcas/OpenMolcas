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

module Trnsfrm

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: IAPR(8), ISP, ISQ, ISR, ISS, LMOP1, LMOQ1, LMOR1, LMOS1, NAP, NAQ, NAR, NAS, NAVX, NBP, NBPQ, NBQ, NBR, NBRS, &
                     NBS, NVXPQ, NX1MX, NX2MX, NX3MX

public :: IAPR, ISP, ISQ, ISR, ISS, LMOP1, LMOQ1, LMOR1, LMOS1, NAP, NAQ, NAR, NAS, NAVX, NBP, NBPQ, NBQ, NBR, NBRS, NBS, NVXPQ, &
          NX1MX, NX2MX, NX3MX

end module Trnsfrm
