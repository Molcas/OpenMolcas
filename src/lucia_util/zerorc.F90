!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1998, Jeppe Olsen                                      *
!***********************************************************************
! Note pt. CC calculations can be restarted from CI
! calculations in the same space by specifying CI=>CC.
! This requires that the input CI vector and the CC
! vector is in the same space. It would be
! better to do the reformatting after the CI.
!
! One can then do f.ex a partial CISDTQ to initialize the
! CCSD.
!
! It would be
!
! Lucia.f : GAS implementing no pair relativistic Theory
!
! Version of March 2000, Jeppe Olsen

subroutine ZERORC(IFIL,IAMPACKED)
! A record was known to be identical zero
!
! Write corresponding info to file IFIL
!
! IAMPACKED added Oct. 98 / Jeppe Olsen

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IFIL, IAMPACKED
integer(kind=iwp) :: ISCR(2)

! Zero record
ISCR(1) = 1
! Packed form
ISCR(2) = IAMPACKED

call ITODS(ISCR,2,2,IFIL)

end subroutine ZERORC
