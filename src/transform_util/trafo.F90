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

module Trafo

! CONTROL DATA FOR TRANSFORMATION SECTION

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: IAD13, IADOUT(3*36*36), ISP, ISQ, ISR, ISS, ITP, ITQ, ITR, ITS, LMOP, LMOP2, LMOQ, LMOQ2, LMOR, LMOR2, LMOS, &
                     LMOS2, LRUPQ, LTUPQ, LURPQ, NBP, NBPQ, NBQ, NBR, NBRS, NBS, NOCP, NOCQ, NOCR, NOCS, NOP, NOQ, NOR, NOS, NPQ

public :: IAD13, IADOUT, ISP, ISQ, ISR, ISS, LMOP, LMOP2, LMOQ, LMOQ2, LMOR, LMOR2, LMOS, LMOS2, LRUPQ, LTUPQ, LURPQ, NBP, NBPQ, &
          NBQ, NBR, NBRS, NBS, NOCP, NOCQ, NOCR, NOCS, NOP, NOQ, NOR, NOS, NPQ
public :: ITP, ITQ, ITR, ITS

end module Trafo
